#!/usr/bin/env python3

"""
This script creates the final DeepOmicsFFPE output.
Author  : DeepOmicsFFPE Team
Date    : 2025-05-13
Version : 1.0.0
Contact : deepomics.ffpe@theragenbio.com
License : © 2025 THERAGEN BIO CO.,LTD. ALL RIGHTS RESERVED.
"""

import argparse
import sys
import os
import pandas as pd
from pathlib import Path
import subprocess
from datetime import datetime
import traceback

def parse_args():
    parser = argparse.ArgumentParser(
        prog='doffpe fin',
        description='\n'
                    'about: Extract features for DeepOmicsFFPE pipeline.\n\n',
        epilog='example: doffpe fin -o output -O DeepOmicsFFPE -p output.pred.tsv',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-ev', '--input-expanded-vcf', required=True, type=Path, help='VCF file with multi-allelic variants split by bcftools norm') # {outdir}/{prefix}.input.expanded.vcf.gz
    parser.add_argument('-o', '--output-prefix', required=True, type=str, help='Prefix to be used for output file names')
    parser.add_argument('-O', '--output-dir', required=False, type=str, default="DeepOmicsFFPE", help='Name of the directory to save the output files')
    parser.add_argument('-p', '--pred-file', required=True, type=Path, help='Inference result file path')
    parser.add_argument('--process-all-variants', action='store_true', help='If specified, include all variants regardless of FILTER status.')

    return parser.parse_args()

def check_files_exist(file_paths):
    missing_files = [path for path in file_paths if not os.path.isfile(path)]
    if missing_files:
        print("\n[ERROR_401] No such file or directory:", file=sys.stderr)
        for file in missing_files:
            print(f" - {file}", file=sys.stderr)
        sys.exit(1)  # exit
    else:
        print("[OK] Ready to start Stage 4: Postprocessing", flush=True)

if __name__ == "__main__":
    
    args = parse_args()
    input_expanded_vcf = args.input_expanded_vcf # f"{outdir}/{prefix}.input.expanded.vcf.gz"
    prefix = args.output_prefix
    outdir = args.output_dir
    pred_tsv = args.pred_file # f"{outdir}/{prefix}.pred.tsv.gz" <- f"{featdir}/{prefix}.pred.tsv.gz"
    process_all_variants = args.process_all_variants

    # START
    stage_start = datetime.now()
    print("\n" + "*"*60)
    print(f"* Stage 4: Postprocessing STARTED at {stage_start.strftime('%Y-%m-%d %H:%M:%S')}")
    print("*"*60 + "\n")
    script_name = os.path.basename(sys.argv[0])
    script_args = sys.argv[1:]
    command_line = ' '.join([script_name] + script_args)
    print(f"* Command: {command_line}\n")
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 0. Chcek the input file
    check_files_exist([input_expanded_vcf]) # File not found error
    
    cmd_fin_vcf = ["bash", f"{script_dir}/ffpe_client/fin_vcf.sh", "-ev", input_expanded_vcf, "-p", pred_tsv, "-o", prefix, "-O", outdir]

    if process_all_variants:
        cmd_fin_vcf = cmd_fin_vcf.append("--process-all-variants")

    try:
        subprocess.run(cmd_fin_vcf, check=True)
        print("\r✅ File saved successfully.")
        
        # Complete
        stage_end = datetime.now()
        print("\n" + "*"*60)
        print(f"* Stage 4: Postprocessing COMPLETED at {stage_end.strftime('%Y-%m-%d %H:%M:%S')}")
        print("*"*60 + "\n")

    except Exception as e:
        print("[ERROR_400] Postprocessing step failed.", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)
