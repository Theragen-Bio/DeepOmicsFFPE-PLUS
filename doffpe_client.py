#!/usr/bin/env python3

"""
DeepOmicsFFPE is used to distinguish somatic variants from formalin-induced artifacts.
Author	: DeepOmicsFFPE Team
Date	: 2025-05-13
Version : 1.0.0
Contact : deepomics.ffpe@theragenbio.com
License : © 2025 THERAGEN BIO CO.,LTD. ALL RIGHTS RESERVED.
"""

import os
import argparse
from datetime import datetime
from pathlib import Path
import subprocess
import sys
import pysam
import traceback
import requests
import json
import copy

# API endpoint URL
base_url = "http://dofp-service.ptbio.kr/api/v1"
analysis_url = f"{base_url}/analysis"
api_key_check_url = f"{base_url}/apikeys/check"
print(os.path.abspath(__file__))

def parse_args():
	parser = argparse.ArgumentParser(
		prog='doffpe',
		description='\n'
					'about: Extract features for DeepOmicsFFPE pipeline.\n'
					'{0}This tool parses VCF files and BAM to produce input-ready features.\n\n'.format(' ' * 8),
		epilog='example: doffpe -v input.vcf -b input.bam -r hg38 -s wgs_pcr -o output -O DeepOmicsFFPE -t 8',
		formatter_class=argparse.RawDescriptionHelpFormatter
	)

	parser.add_argument('-v', '--variant-file', required=True, type=Path, help='Input variants file path')
	parser.add_argument('-b', '--bam-file', required=True, type=Path, help='Input BAM file path')
	parser.add_argument('-r', '--ref-version', required=True, choices=['hg19', 'hg38'], help='Reference genome version: hg19 or hg38')
	parser.add_argument('-s', '--seq-type', required=True, choices=['wes', 'wgs_pcr', 'wgs_pcrfree'], help='Sequencing type: "wes" for Whole Exome Sequencing, "wgs_pcr" for Whole Genome Sequencing with PCR-based library prep, or "wgs_pcrfree" for Whole Genome Sequencing with PCR-free library prep.')
	parser.add_argument('-o', '--output-prefix', required=True, type=str, help='Prefix to be used for output file names')
	parser.add_argument('-O', '--output-dir', required=False, type=str, default="DeepOmicsFFPE", help='Name of the directory to save the output files')
	parser.add_argument('-t', '--threads', required=False, type=int, default=0, help='Use multithreading with <int> worker threads')
	parser.add_argument('--process-all-variants', action='store_true', help='If specified, include all variants regardless of FILTER status.')
	parser.add_argument('--api-key', required=False, default = None, help="To use this program, you must provide an API token. If you have used it previously, the token may already be stored in the [home_directory]/.dofp path, and you won't need to provide it again.")
	parser.add_argument('--environment', default = None, help='If not specified, the script will automatically be executed in the dofp environment.')
	
	return parser.parse_args()

def check_files_exist(file_paths):
	missing_files = [path for path in file_paths if not os.path.isfile(path)]
	if missing_files:
		print("\n[ERROR_001] No such file or directory:", file=sys.stderr)
		for file in missing_files:
			print(f" - {file}", file=sys.stderr)
		sys.exit(1)  # exit
	else:
		print("[OK] All files exist.", flush=True)


def validate_api_key(api_key):
	headers = {
		'Authorization': 'ApiKey {0}'.format(api_key)
	}

	# 1. set home directory
	home_dir = os.path.expanduser("~")
	dofp_dir = os.path.join(home_dir, ".dofp")
	credentials_path = os.path.join(dofp_dir, "credentials")

	credential = {}
	if api_key == None :
		if os.path.isfile(credentials_path) == True : 
			with open(credentials_path, "r", encoding="utf-8") as f:
				credential = json.load(f)
		else : 
			print("Please use --api-key option to initialize API key. Once key is enrolled, You don't need to use --api-key option.")
			sys.exit()
		# Check API key
		if credential != {} and "api_key" in credential : 
			api_key = credential["api_key"]
		else : 
			print("Credential file is corrupted. Please use again --api-key option.")
			sys.exit()
	else : 
		response = requests.get(api_key_check_url, headers=headers)
		print("validate_api_key response.text", response.text)
		data = response.json()
		print("data", data)
		if data['isValid'] == False : 
			print(data['message'])
			sys.exit()

	# 2. create ~/.dofp folder
	os.makedirs(dofp_dir, exist_ok=True)

	# 3. create credentials file (json)
	credentials = {"api_key": api_key}
	with open(credentials_path, "w") as f:
		json.dump(credentials, f, indent=4)
		print(f"✅ API key saved to {credentials_path}")
	
	return api_key


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
	api_key = args.api_key
	
	if outdir is None:
		outdir = os.path.abspath(os.getcwd()) + "/"
	else:
		outdir = (outdir + "/").replace("//", "/")

	api_key = validate_api_key(api_key)
	
	# START
	stage_start = datetime.now()
	print("\n" + "*"*60)
	print(f"* DeepOmicsFFPE STARTED at {stage_start.strftime('%Y-%m-%d %H:%M:%S')}")
	print("*"*60 + "\n")
	script_name = os.path.basename(sys.argv[0])
	script_args = sys.argv[1:]
	command_line = ' '.join([script_name] + script_args)
	print(f"* Command: {command_line}\n")

	# 1. Chcek the input file
	check_files_exist([input_vcf]) # File not found error
	
	## check the bam index file
	with pysam.AlignmentFile(input_bam, "rb") as bamfile:
		if not bamfile.has_index():
			print("\n[ERROR_002] BAM Index file does not exist.", file=sys.stderr)
			sys.exit(1)
		else:
			print("[OK] BAM Index file exists.", flush=True)
	
	## DeepOmicsFFPE file path
	input_expanded_vcf = f"{outdir}{prefix}.input.expanded.vcf.gz"
	allele_file = f"{outdir}{prefix}.allele_data.pkl.gz"
	pred_tsv = f"{outdir}test_result.tsv.gz"
	
	# 2. Run
	base_dir = os.path.dirname(os.path.abspath(__file__))
	script_dir = os.path.join(base_dir, "ffpe_utils")
	
	cmd_01 = ["python", f"{script_dir}/01_ext_bam.py", "-v", input_vcf, "-b", input_bam, "-r", ref_ver, "-s", seq_type, "-o", prefix, "-O", outdir, "-t", str(num_proc)]
	cmd_04 = ["python", f"{script_dir}/04_fin.py", "-ev", input_expanded_vcf, "-o", prefix, "-O", outdir, "-p", pred_tsv]

	if process_all_variants:
		cmd_01 = cmd_01.append("--process-all-variants")
		cmd_04 = cmd_04.append("--process-all-variants")

	try:
		subprocess.run(cmd_01, check=True)
	except Exception as e:
		print(f"\n[ERROR_100] Data PreProcessing step failed.\n", file=sys.stderr)
		traceback.print_exc()
		sys.exit(1) # error

	headers = {
		'Authorization': 'ApiKey {0}'.format(api_key)
	}

	# JSON parameter
	parameters = {
		"variant-file": f"{input_vcf}",
		"bam-file": f"{input_bam}",
		"ref-version": f"{ref_ver}",
		"seq-type": f"{seq_type}",
		"output-prefix": f"{prefix}",
		"output-dir": f"{outdir}",
		"threads": num_proc,
		"process-all-variants": f"{process_all_variants}", 
		"variant-read-counts": os.path.basename(allele_file)
	}

	files = {
		'file': (os.path.basename(allele_file), open(f'{allele_file}', 'rb'), 'application/octet-stream')
	}

	data = {
		'script': 'doffpe_server',
		'language': 'python',
		'envType': 'conda',
		'envPath': '/dofp-data/env/conda/miniconda3',
		'envName': 'dofp',
		'parameters': json.dumps(parameters)
	}

	# request analysis api
	response = requests.post(analysis_url, headers=headers, files=files, data=data)

	print(response.status_code)
	print("response.text", response.text)
	data = response.json()
	print("data", data)

	# log stream
	url_stream = f"{analysis_url}/{data['seq']}/logs/stream"
	headers_stream = copy.deepcopy(headers)
	headers_stream['Accept'] = "text/event-stream"

	# stream=True for log streaming
	resp = requests.get(url_stream, headers=headers_stream, stream=True)
	resp.raise_for_status()

	for raw_line in resp.iter_lines(decode_unicode=False):
		if raw_line:  # To skip empty lines
			try : 
				line = json.loads(raw_line.decode('utf-8', errors='replace')[len('data:'):])['data']
				print(f"date : {line['createdDate']}")
				print(f"analysis step : {line['analysisStep']}")
				print(f"log : {line['log'].strip()}")
				print(f"status : {line['status']}")
				print("")
			except Exception as e : 
				pass
				#line = raw_line.decode('utf-8', errors='replace')
				#print("Received:", line)

	resp.close()

	url_download = f"{analysis_url}/{data['seq']}/result/download"

	# GET intermediate output
	resp = requests.get(url_download, headers=headers)
	resp.raise_for_status()

	# save intermediate output
	with open(pred_tsv, "wb") as f:
		f.write(resp.content)

	print(f"Complete file download : {pred_tsv}")

	# post analysis for client-side
	try:
		subprocess.run(cmd_04, check=True)
	except Exception as e:
		print(f"\n[ERROR_400] Postprocessing step failed.\n", file=sys.stderr)
		traceback.print_exc()
		sys.exit(1) # error
	   
	# Complete
	stage_end = datetime.now()
	print("\n" + "*"*60)
	print(f"* DeepOmicsFFPE COMPLETED at {stage_end.strftime('%Y-%m-%d %H:%M:%S')}")
	print("*"*60 + "\n")
	sys.stdout.flush()
	# sys.stderr.flush()
