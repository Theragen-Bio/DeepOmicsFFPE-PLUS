#!/bin/bash

# ---------------------------------------------------------------------------
# This script creates the final DeepOmicsFFPE output
# for variant files in VCF or compressed VCF (.vcf.gz) format.
# Author  : DeepOmicsFFPE Team
# Date    : 2025-05-13
# Version : 1.0.0
# Contact : deepomics.ffpe@theragenbio.com
# License : © 2025 THERAGEN BIO CO.,LTD. ALL RIGHTS RESERVED.
# ---------------------------------------------------------------------------

# Initialize variables
print_help() {
    echo
    echo Usage: bash fin_vcf.sh [OPTIONS]
    echo
    echo Options:
    echo  "-ev, --input-expanded-vcf"
    echo  "-p, --pred-file"
	echo  "-o, --output-prefix"
    echo  "-O, --output-dir"
    echo  "--process-all-variants"
    echo  -h, --help Display this help message and exit.
    echo
}
echo About to parse command line arguments...


# Parse command line arguments
process_all_variants=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -ev|--input-expanded-vcf) input_expanded_vcf=$2; shift ;;
		-p|--pred-file) pred_tsv=$2; shift ;;
		-o|--output-prefix) prefix=$2; shift ;;
        -O|--output-dir) outdir=$2; shift ;;
        --process-all-variants) process_all_variants=true ;;
        -h|--help) print_help; exit 0 ;;
        *) echo Unknown parameter passed: $1; print_help; exit 1 ;;
    esac
    shift
done
echo Finished parsing arguments.

# 0. Check required arguments
if [[ -z "$input_expanded_vcf" || -z "$pred_tsv" || -z "$prefix" || -z "$outdir" ]]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 -ev input.expanded.vcf.gz -o output -O DeepOmicsFFPE -p output.DeepOmicsFFPE.tsv"
    exit 1
fi
echo All parameters seem okay. About to run DeepOmicsFFPE command...

# 1. Ready to files
vcf_header=$'##fileformat=VCFv4.2\n'\
$'##DeepOmicsFFPE_Version=0.0.1\n'\
$'##INFO=<ID=DeepOmicsFFPE_score,Number=1,Type=Float,Description="Predicted FFPE artifact score">\n'\
$'##INFO=<ID=IS_VARIANT,Number=1,Type=Integer,Description="Predicted FFPE artifact label (0 = Artifact or 1 = True variant)">\n'\
$'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'

## 1-1. DeepOmicsFFPE results (pred.tmp.vcf.gz)
pred_vcf_tmp=${pred_tsv/.tsv.gz/.tmp.vcf.gz}

score_col=$(zcat ${pred_tsv} | head -1 | tr '\t' '\n' | grep -n '^DeepOmicsFFPE_score$' | cut -d: -f1)
pred_col=$(zcat ${pred_tsv} | head -1 | tr '\t' '\n' | grep -n '^IS_VARIANT$' | cut -d: -f1)

{
  echo "$vcf_header"
  zcat "${pred_tsv}" | awk -v s=${score_col} -v p=${pred_col} '
    NR > 1 {
      score = int($s * 1000) / 1000
      printf "%s\t%s\t.\t%s\t%s\t.\t.\tDeepOmicsFFPE_score=%.3f;IS_VARIANT=%d\n", $1, $2, $3, $4, score, $p
    }
  '
} | bgzip -c > "${pred_vcf_tmp}" && tabix -p vcf "${pred_vcf_tmp}"

## 1-2. merge pass & non-pass variants (pred_vcf_tmp & input_expanded_vcf_npass_tmp)
pred_vcf_merge_tmp=${pred_tsv/.tsv.gz/.tmp.merge.vcf.gz}

if [ "$process_all_variants" = true ]; then
    mv ${pred_vcf_tmp} ${pred_vcf_merge_tmp}
else
    input_expanded_vcf_npass_tmp=${input_expanded_vcf/.vcf.gz/.npass.vcf.gz}
    {
    echo "$vcf_header"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -i 'FILTER!="PASS"' "${input_expanded_vcf}" \
        | awk 'BEGIN{OFS="\t"} {print $1, $2, ".", $3, $4, ".", ".", "DeepOmicsFFPE_score=.;IS_VARIANT=."}'
    } | bgzip -c > ${input_expanded_vcf_npass_tmp} && tabix -p vcf ${input_expanded_vcf_npass_tmp}

    bcftools merge ${pred_vcf_tmp} ${input_expanded_vcf_npass_tmp} -Oz -o ${pred_vcf_merge_tmp} && tabix -p vcf ${pred_vcf_merge_tmp}
fi

# 2. Annotate the DeepOmicsFFPE results
doffpe_vcf=${outdir}/${prefix}.DeepOmicsFFPE.vcf.gz
tabix -p vcf ${input_expanded_vcf}
bcftools annotate -a ${pred_vcf_merge_tmp} -c CHROM,POS,REF,ALT,.INFO/DeepOmicsFFPE_score,.INFO/IS_VARIANT -o ${doffpe_vcf} -O z --no-version ${input_expanded_vcf}

# 3. Save the results
doffpe_true_vcf=${outdir}/${prefix}.DeepOmicsFFPE.filtered.vcf.gz
bcftools view -i 'INFO/IS_VARIANT=1' ${doffpe_vcf} -o ${doffpe_true_vcf} -O z

# 4. Remove the temp files
rm ${pred_vcf_tmp} ${input_expanded_vcf_npass_tmp} ${pred_vcf_merge_tmp} ${pred_vcf_tmp}.tbi ${input_expanded_vcf_npass_tmp}.tbi ${pred_vcf_merge_tmp}.tbi
rm ${input_expanded_vcf}.tbi