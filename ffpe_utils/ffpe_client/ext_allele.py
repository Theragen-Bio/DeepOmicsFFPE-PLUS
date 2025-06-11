#!/usr/bin/env python3

import re

class Variant:
    max_read_length = 0
    mapping_quality_threshold = 20
    def __init__(
            self, 
            bam_object, 
            chrom: str, 
            var_pos: int, 
            ref: str, 
            alt: str, 
            ):
        """
        parameters
        bam_object : pysam.AlignmentFile object
        chrom : chromosome
        var_pos : position of variant
        ref : reference sequence
        alt : variant sequence
        """
        self.bam = bam_object
        self.chrom = chrom
        self.var_pos = var_pos
        self.ref = ref
        self.alt = alt
        self.pseudo_count = 1

        self.dp = 0
        self.ad_ref = 0
        self.ad_alt= 0
        self.ad_other= 0

        self.ref_sup_r1 = 0
        self.ref_sup_r2 = 0
        self.alt_sup_r1 = 0
        self.alt_sup_r2 = 0

        self.f1r2_ref = 0
        self.f2r1_ref = 0
        self.f1r2_alt = 0
        self.f2r1_alt = 0

        self.sb_1_ref = 0
        self.sb_2_ref = 0
        self.sb_1_alt = 0
        self.sb_2_alt = 0

        self.isize_ref_list = [] # list of length of insert supporting reference allele
        self.isize_alt_list= [] 
        self.map_error_rate_ref_list = []
        self.map_error_rate_alt_list = []
        self.map_qualities_ref_list = []
        self.map_qualities_alt_list = []
        self.dist_to_end_ref_list = []
        self.dist_to_end_alt_list = []

        self.ref_support_reads = []
        self.alt_support_reads = []
        self.other_support_reads = []

        self.ref_qname = []
        self.alt_qname = []
        self.other_qname = []
        self.ambigous_qname = []

        self._get_metrics()

    def _get_genomic_coordinates(self, read):
        expanded_cigarstring = [char for count, char in re.findall(r'(\d+)(\D)', read.cigarstring) for char in int(count) * char]
        softclip_idx = [i for i, char in enumerate(expanded_cigarstring) if char == "S"]
        genomic_coordinates_rm_softclip = [pos for i, pos in enumerate(read.get_reference_positions(full_length=True)) if i not in softclip_idx]
        genomic_coordinates_1based = [pos + 1 if pos is not None else None for pos in genomic_coordinates_rm_softclip]
        return genomic_coordinates_1based



    def _ref_alt_other_support(self, read):
        """ 
        return values:
        0 means that the read supports ref allele
        1 means that the read supports alt allele
        2 means that the read does not support neither ref nor alt
        """
        seq = read.query_alignment_sequence
        genomic_coordinates = self._get_genomic_coordinates(read)
        # genomic_coordinates is a list of genomic position of every base in the read. 
        # genomic_coordinates can include None values corresponding to insertions

        if self.var_pos not in genomic_coordinates:
            return 2
        else:
            idx_for_var_pos = genomic_coordinates.index(self.var_pos)

        if len(self.ref) == len(self.alt): # SNV
            """
            In case of SNV, 
            chr1    10  A  G
            pos_of_interest is going to be list(range(10, 11)) -> [10]
            """
            pos_of_interest = list(range(self.var_pos, self.var_pos + len(self.ref)))
            if any(p not in genomic_coordinates for p in pos_of_interest):
                # any or all ?????????????????????? 
                return 2
            else:
                if seq[idx_for_var_pos : idx_for_var_pos + len(self.alt)] == self.alt:
                    return 1
                elif seq[idx_for_var_pos : idx_for_var_pos + len(self.ref)] == self.ref:
                    return 0
                else:
                    return 2

        elif len(self.ref) > len(self.alt): # Deletion
            """
            In case of deletion, 
            chr1    10  ATT A
            pos_of_interest is going to be [11, 12]
            """
            pos_of_interest = list(range(self.var_pos + 1, self.var_pos + len(self.ref)))
            if pos_of_interest[-1] + 1 not in genomic_coordinates:
                return 2
            else:
                if all(p not in genomic_coordinates for p in pos_of_interest): 
                    if seq[idx_for_var_pos : idx_for_var_pos + len(self.alt)] == self.alt:
                        return 1
                    else:
                        return 2
                elif all(p in genomic_coordinates for p in pos_of_interest):
                    if seq[idx_for_var_pos : idx_for_var_pos + len(self.ref)] == self.ref:
                        return 0
                    else:
                        return 2
                else:
                    return 2

        else: # Insertion len(self.ref) < len(self.alt)
            len_insert = len(self.alt) - len(self.ref)
            if (self.var_pos + 1) not in genomic_coordinates:
                return 2
            else:
                if (self.var_pos + 1) == genomic_coordinates[idx_for_var_pos + 1]:
                    if seq[idx_for_var_pos : idx_for_var_pos + len(self.ref)] == self.ref:
                        return 0
                    else:
                        return 2
                elif [None] * len_insert == genomic_coordinates[idx_for_var_pos + 1 : idx_for_var_pos + 1 + len_insert]:
                    if seq[idx_for_var_pos : idx_for_var_pos + len(self.alt)] == self.alt:
                        return 1
                    else:
                        return 2
                else:
                    return 2
    
    def _is_valid(self, read):
        if read.mapping_quality > Variant.mapping_quality_threshold and \
                read.is_proper_pair and \
                (not read.is_qcfail) and \
                (not read.is_duplicate) and \
                (not read.is_secondary) and \
                (not read.is_supplementary) and \
                self.var_pos - 1 in read.positions:
            return True
        else:
            return False

    def _get_dist_to_end(self, read):
        # if read.is_forward: # iykim edit 25.04.07
        if not read.is_reverse:
            dist = read.positions.index(self.var_pos-1)
        else:
            dist = len(read.positions) - read.positions.index(self.var_pos-1) - 1
        return dist

    def _update_information(self):
        for record in self.ref_support_reads:
            if record['name'] in self.ambigous_qname: continue
            self.f1r2_ref += record['f1r2']
            self.f2r1_ref += record['f2r1']
            self.sb_1_ref += record['sb_1']
            self.sb_2_ref += record['sb_2']
            self.ref_sup_r1 += record['sup_r1']
            self.ref_sup_r2 += record['sup_r2']
            self.map_qualities_ref_list.append(record['mapping_quality'])
            self.dist_to_end_ref_list.append(record['dist_to_end'])
            try:
                self.isize_ref_list.append(record['isize'])
            except KeyError:
                pass
        for record in self.alt_support_reads:
            if record['name'] in self.ambigous_qname: continue
            self.f1r2_alt += record['f1r2']
            self.f2r1_alt += record['f2r1']
            self.sb_1_alt += record['sb_1']
            self.sb_2_alt += record['sb_2']
            self.alt_sup_r1 += record['sup_r1']
            self.alt_sup_r2 += record['sup_r2']
            self.map_qualities_alt_list.append(record['mapping_quality'])
            self.dist_to_end_alt_list.append(record['dist_to_end'])
            try:
                self.isize_alt_list.append(record['isize'])
            except KeyError:
                pass


    def to_dict(self):
        self._update_information()

        tmp = dict()
        tmp['chrom'] = self.chrom
        tmp['pos'] = self.var_pos
        tmp['ref'] = self.ref
        tmp['alt'] = self.alt 
        tmp['max_read_length'] = Variant.max_read_length
        tmp['dp'] = self.dp
        tmp['ad_ref'] = len(set([name for name in self.ref_qname if name not in self.ambigous_qname]))
        tmp['ad_alt'] = len(set([name for name in self.alt_qname if name not in self.ambigous_qname]))
        tmp['ad_other'] = len(set([name for name in self.other_qname if name not in self.ambigous_qname]))
        tmp['ref_sup_r1'] = self.ref_sup_r1
        tmp['ref_sup_r2'] = self.ref_sup_r2
        tmp['alt_sup_r1'] = self.alt_sup_r1
        tmp['alt_sup_r2'] = self.alt_sup_r2
        tmp['f1r2_ref'] = self.f1r2_ref
        tmp['f2r1_ref'] = self.f2r1_ref
        tmp['f1r2_alt'] = self.f1r2_alt
        tmp['f2r1_alt'] = self.f2r1_alt
        tmp['sb_1_ref'] = self.sb_1_ref
        tmp['sb_2_ref'] = self.sb_2_ref
        tmp['sb_1_alt'] = self.sb_1_alt
        tmp['sb_2_alt'] = self.sb_2_alt
        tmp['isize_ref'] = self.isize_ref_list
        tmp['isize_alt'] = self.isize_alt_list
        tmp['map_qualities_ref'] = self.map_qualities_ref_list
        tmp['map_qualities_alt'] = self.map_qualities_alt_list
        tmp['dist_to_end_ref'] = self.dist_to_end_ref_list
        tmp['dist_to_end_alt'] = self.dist_to_end_alt_list
        return tmp



    def _get_metrics(self):
        for read in self.bam.fetch(self.chrom, self.var_pos-1, self.var_pos):
            if not self._is_valid(read):
                continue
            else:
                supporting_allele = self._ref_alt_other_support(read)
                if Variant.max_read_length < len(read.query_alignment_sequence):
                    Variant.max_read_length = len(read.query_alignment_sequence)
            
            if supporting_allele == 0: # support ref allele
                if read.qname in self.alt_qname or read.qname in self.other_qname:
                    self.ambigous_qname.append(read.qname) # inconsistently support allele
                    self.dp -= 1
                    continue
                if read.qname not in self.ref_qname:
                    self.dp += 1
                    read_dict = {'name': read.qname, 'allele': 'ref', 'isize': read.isize}
                    self.ref_qname.append(read.qname)
                else:
                    read_dict = {'name': read.qname, 'allele': 'ref'}

                # if (read.is_forward and read.is_read1) or (read.is_reverse and read.is_read2): # F1R2 # iykim edit 25.04.07
                if (not read.is_reverse and read.is_read1) or (read.is_reverse and read.is_read2): # F1R2
                    read_dict['f1r2'] = 1
                    read_dict['f2r1'] = 0
                else:
                    read_dict['f1r2'] = 0
                    read_dict['f2r1'] = 1
                # if read.is_forward: # SB1 # iykim edit 25.04.07
                if not read.is_reverse: # SB1
                    read_dict['sb_1'] = 1
                    read_dict['sb_2'] = 0
                else:
                    read_dict['sb_1'] = 0
                    read_dict['sb_2'] = 1
                if read.is_read1:
                    read_dict['sup_r1'] = 1
                    read_dict['sup_r2'] = 0
                else:
                    read_dict['sup_r1'] = 0
                    read_dict['sup_r2'] = 1
                read_dict['mapping_quality'] = read.mapping_quality
                read_dict['dist_to_end'] = self._get_dist_to_end(read)
                self.ref_support_reads.append(read_dict)


            elif supporting_allele == 1: # support alt allele
                if read.qname in self.ref_qname or read.qname in self.other_qname:
                    self.ambigous_qname.append(read.qname) # inconsistently support allele
                    self.dp -= 1
                    continue
                if read.qname not in self.alt_qname:
                    self.dp += 1
                    read_dict = {'name': read.qname, 'allele': 'alt', 'isize': read.isize}
                    self.alt_qname.append(read.qname)
                else:
                    read_dict = {'name': read.qname, 'allele': 'alt'}

                # if (read.is_forward and read.is_read1) or (read.is_reverse and read.is_read2): # F1R2 # iykim edit 25.04.07
                if (not read.is_reverse and read.is_read1) or (read.is_reverse and read.is_read2): # F1R2
                    read_dict['f1r2'] = 1
                    read_dict['f2r1'] = 0
                else:
                    read_dict['f1r2'] = 0
                    read_dict['f2r1'] = 1
                # if read.is_forward: # SB1 # iykim edit 25.04.07
                if not read.is_reverse: # SB1
                    read_dict['sb_1'] = 1
                    read_dict['sb_2'] = 0
                else:
                    read_dict['sb_1'] = 0
                    read_dict['sb_2'] = 1
                if read.is_read1:
                    read_dict['sup_r1'] = 1
                    read_dict['sup_r2'] = 0
                else:
                    read_dict['sup_r1'] = 0
                    read_dict['sup_r2'] = 1
                read_dict['mapping_quality'] = read.mapping_quality
                read_dict['dist_to_end'] = self._get_dist_to_end(read)
                self.alt_support_reads.append(read_dict)
            else:
                if read.qname in self.ref_qname or read.qname in self.alt_qname:
                    self.ambigous_qname.append(read.qname) # inconsistently support allele
                    self.dp -= 1
                    continue
                if read.qname not in self.other_qname:
                    self.dp += 1
                    read_dict = {'name': read.qname, 'allele': 'other'}
                    self.other_qname.append(read.qname)
                self.other_support_reads.append(read_dict)
