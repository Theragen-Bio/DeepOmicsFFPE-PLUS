#!/usr/bin/env python3

"""
This script extracts information from a BAM file for use in DeepOmicsFFPE analysis.
Author  : DeepOmicsFFPE Team
Date    : 2025-05-13
Version : 1.0.0
Contact : deepomics.ffpe@theragenbio.com
License : © 2025 THERAGEN BIO CO.,LTD. ALL RIGHTS RESERVED.
"""

import os
import time
import argparse
from pathlib import Path
import concurrent.futures
from itertools import chain
import pandas as pd
from math import ceil
import pysam
import io
import gzip
from tqdm import tqdm
from ffpe_client import ext_allele
import threading
import pickle
import sys
from datetime import datetime
import traceback
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(
        prog='ext_bam',
        description='\n'
                    'about: Extract features for DeepOmicsFFPE pipeline.\n'
                    '       This tool parses VCF/TSV/CSV files and BAM to produce input-ready features.\n\n',
        epilog='example: ext_bam -v input.vcf -b input.bam -r hg38 -s wgs_pcr -o output -O DeepOmicsFFPE -t 8',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-v', '--variant-file', required=True, type=Path, help='Input variants file path')
    parser.add_argument('-b', '--bam-file', required=True, type=Path, help='Input BAM file path')
    parser.add_argument('-r', '--ref-version', required=True, choices=['hg19', 'hg38'], help='Reference genome version: hg19 or hg38')
    parser.add_argument('-s', '--seq-type', required=True, choices=['wes', 'wgs_pcr', 'wgs_pcrfree'], help='Sequencing type: "wes" for Whole Exome Sequencing, "wgs_pcr" for Whole Genome Sequencing with PCR-based library prep, or "wgs_pcrfree" for Whole Genome Sequencing with PCR-free library prep.')
    parser.add_argument('-o', '--output-prefix', required=True, type=str, help='Prefix to be used for output file names')
    parser.add_argument('-O', '--output-dir', required=False, type=str, default='DeepOmicsFFPE', help='Name of the directory to save the output files')
    parser.add_argument('-t', '--threads', required=False, type=int, default=0, help='Use multithreading with <int> worker threads')

    # parser.add_argument('--predict-only-passed-calls', required=True, choices=['True', "False"], default='True', help='If True, process only PASS variants; otherwise, process all.')
    parser.add_argument('--process-all-variants', action='store_true', help='If specified, include all variants regardless of FILTER status.')

    return parser.parse_args()

def read_vcf_gz(path):
    with gzip.open(path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')]
    df = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
    df.columns = df.columns.str.lower()
    return df
    
def read_file(path, input_type='vcf.gz'):
    df = read_vcf_gz(path)
    ## check the columns error
    required_cols = ['chrom', 'pos', 'ref', 'alt', 'filter']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"\n[ERROR_102] The variant file is missing required columns.", file=sys.stderr)
        print(f"... Variant file: {path}", file=sys.stderr)
        print(f"... Missing columns: {missing_cols}", file=sys.stderr)
        sys.exit(1)  # error
    else:
        print("[OK] All required columns are present in the variant file.", flush=True)
    
    return df[required_cols]

def extract_allele_depth(df_chunk, bam_path):
    bam = pysam.AlignmentFile(bam_path)
    results = []
    for _, row in df_chunk.iterrows():
        variant = ext_allele.Variant(bam, row['chrom'], row['pos'], row['ref'], row['alt'])
        results.append(variant.to_dict())
    return results

def extract_allele_depth_wrapper(args): # multi-process optimization
        return extract_allele_depth(*args)

def show_saving_spinner(stop_event):
    spinner = ['|', '/', '-', '\\']
    idx = 0
    while not stop_event.is_set():
        print(f"\rAnalyzing... {spinner[idx % len(spinner)]}", end="", flush=True)
        idx += 1
        time.sleep(0.1)
    print('\rAnalyzing... Done!     ', flush=True)

if __name__ == "__main__":

    args = parse_args()
    input_vcf = args.variant_file
    input_bam = args.bam_file
    ref_ver = args.ref_version
    seq_type = args.seq_type
    prefix = args.output_prefix
    outdir = args.output_dir
    num_proc = args.threads
    process_all_variants = args.process_all_variants

    # START
    stage_start = datetime.now()
    print("\n" + "*"*60)
    print(f"* Stage 1: Data Preprocessing STARTED at {stage_start.strftime('%Y-%m-%d %H:%M:%S')}")
    print("*"*60 + "\n")
    script_name = os.path.basename(sys.argv[0])
    script_args = sys.argv[1:]
    command_line = ' '.join([script_name] + script_args)
    print(f"* Command: {command_line}\n")

    # 1. Create a directory to store the results
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    ## 1-3. Expanding multi-allelic variants
    print(f"\n✅ Expanding multi-allelic variants")
    try:
        input_expanded_vcf = f"{outdir}/{prefix}.input.expanded.vcf.gz"
        cmd_multi_allele = ["bcftools", "norm", "-Ov", "-m-any", "--force", input_vcf, "--output", input_expanded_vcf]
        subprocess.run(cmd_multi_allele, check=True)
    except Exception as e:
        print("\n[ERROR_101] Please check the 'variant file path' & 'VCF file format' again.:", file=sys.stderr)
        print(f"... File name: {input_vcf}", file=sys.stderr)
        print("The file must include the following columns: [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO]", file=sys.stderr)
        print("Supported file formats are: [.vcf, .vcf.gz]\n", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1) # error

    # 2. Load the input VCF file.
    print(f"\n✅ Load the variants file: {input_expanded_vcf}")
    try:
        vcf_df = read_file(input_expanded_vcf)
        print(f"\r✅ Number of Variants: {len(vcf_df):,}")
    except Exception as e:
        print("\n[ERROR_102] Please check the 'variant file path' & 'VCF file format' again.:", file=sys.stderr)
        print(f"... File name: {input_expanded_vcf}", file=sys.stderr)
        print("The file must include the following columns: [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO]", file=sys.stderr)
        print("Supported file formats are: [.vcf, .vcf.gz]\n", file=sys.stderr)
        sys.exit(1) # variant file type error
    
    ## 2-1. Check the chrom format: ex. chr1 ...
    chroms = vcf_df['chrom'].unique()
    if not all(str(chrom).startswith('chr') for chrom in chroms):
        raise ValueError(
            "\n[ERROR_103] The 'CHROM' column contains invalid values.\n"
            f"... File name: {input_vcf}\n"
            "Expected chromosome names to start with 'chr' (e.g., 'chr1', 'chr2', ..., 'chrX', 'chrY').\n"
            "Please verify and correct the input file format.\n"
        )
    
    ## 2-2. FILTER variants: process_all_variants
    if process_all_variants:
        print(f"Processing all {len(vcf_df)} variants (including non-PASS).")
    else:
        ## 2-3. Check the FILTER: All values in the 'FILTER' column are '.' -> error
        if vcf_df['filter'].eq('.').all():
            raise ValueError(
                "\n[ERROR_104] All values in the 'FILTER' column are '.'.\n"
                f"... File name: {input_vcf}\n"
                "\nThis error indicates that the VCF file may not have been filtered, or the filtering step did not populate the FILTER column properly\n."
                "\n1. We recommend checking the filtering settings of your variant caller, or applying a post-processing filter to the VCF file.\n"
                "Our pipeline is designed to analyze only variants with `FILTER=PASS`, and we can only guarantee results based on such filtered variants.\n"
                "\n2. If you prefer to analyze all variants regardless of the FILTER status, you may use the `--process-all-variants` option.\n"
                "However, please note that we do not recommend analyzing variants that have not passed standard quality filters, as this may affect data reliability.\n"
            )
        
        # vcf_df_nonpass = vcf_df[vcf_df['filter'] != 'PASS'].copy()
        # vcf_df['to_analyze'] = (vcf_df['FILTER'] == 'PASS').astype(int)
        vcf_df = vcf_df[vcf_df['filter'] == 'PASS'].copy()
        print(f"Processing {len(vcf_df)} PASS variants only.")

    # 3. Manage num of workers and split dataframe accordingly
    ## workers
    max_num_worker = os.cpu_count()  # All of CPU

    if num_proc == 0 or num_proc >= max_num_worker:
        num_worker = os.cpu_count()  # All of CPU
    else:
        num_worker = num_proc
    num_chunks = num_worker * 3

    ## dataframe
    vcf_df = vcf_df.reset_index(drop=True)
    bs = int(ceil(len(vcf_df) / num_chunks))
    idx_for_split_df = vcf_df.index.to_list()[::bs]
    idx_for_split_df.append(vcf_df.index[-1] + 1)

    idx_pairs = [(idx_for_split_df[i], idx_for_split_df[i+1]) for i in range(len(idx_for_split_df) -1)]  # index pairs: (start, end)
    args_list = [(vcf_df[start:end], input_bam) for start, end in idx_pairs]  # Generate a list of (VCF chunk, BAM file) pairs based on index ranges # input of multi-process

    # 4. Extract raw read data from bam # multi-processing
    start_time = time.perf_counter()

    ## Display processing status in the CLI
	#stop_event = threading.Event()
	#spinner_thread = threading.Thread(target=show_saving_spinner, args=(stop_event,))
	#spinner_thread.start()

    try:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_worker) as executor:
            result = executor.map(extract_allele_depth_wrapper, args_list)
            allele_data = list(chain(*tqdm(result, total=len(args_list), desc="✅ Extracting data from BAM file")))

        # with concurrent.futures.ProcessPoolExecutor(max_workers=num_worker) as executor:
        #     result = executor.map(extract_allele_depth_wrapper, args_list)
        #     result = list(result)  # 여기서 모든 작업이 끝날 때까지 대기
        #     allele_data = list(chain(*tqdm(result, total=len(args_list), desc="✅ Extracting data from BAM file")))

        
        ## Check the results
        print(f"\n✅ Total allele records: {len(allele_data):,}")
        vnum_bef = len(vcf_df)
        vnum_aft = len(allele_data)

        if vnum_bef != vnum_aft:
            sys.stderr.write(f"[ERROR_105] Data PreProcessing step failed. Variant count mismatch: before={vnum_bef}, after={vnum_aft}\n")
            sys.exit(1)
        
        after_bam_time = time.perf_counter()
        
        # 5. Save the allele data
        print("\r✅ Save the allele data")
        save_start_time = time.perf_counter()

        ## Save the pickle file
        with gzip.open(f"{outdir}/{prefix}.allele_data.pkl.gz", "wb") as f:
            pickle.dump(allele_data, f, protocol=pickle.HIGHEST_PROTOCOL)
        
        end_time = time.perf_counter()

    except Exception as e:
        print("\n[ERROR_100] Data PreProcessing step failed.\n", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)

    finally:
		#stop_event.set()  # Signal that the task is complete
		#spinner_thread.join()  # Wait until the spinner thread finishes
    
        print("\r✅ File saved successfully.")

        # Complete
        stage_end = datetime.now()
        print("\n" + "*"*60)
        print(f"* Stage 1: Data Preprocessing COMPLETED at {stage_end.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"\r✅ Number of cpu: {num_worker}, Sample: {prefix}, Number of Variants: {len(vcf_df):,}, For extracting info from bam: {round(after_bam_time - start_time, 2)}s, For Saving file: {round(end_time - save_start_time, 2)}s")
        print("*"*60 + "\n")
