#!/usr/bin/env python3
import os
import requests
import argparse
import time
import sys
import json
import subprocess
from requests_toolbelt.multipart.encoder import MultipartEncoder

def run_alphafold_local(fasta_path, output_dir, db_path):
    cmd = [
        "python3 /app/alphafold/run_alphafold.py",
        f"--fasta_paths={fasta_path}",
        f"--output_dir={output_dir}",
        f"--data_dir={db_path}",
        "--max_template_date=2024-01-01",
        "--model_preset=monomer",
        "--db_preset=full_dbs"
    ]
    subprocess.run(" ".join(cmd), shell=True, check=True)

def run_alphafold(msa_file, id, local_db, fasta_file):
    if local_db and os.path.exists(local_db):
        print(f"Local DB found at {local_db}. Running local AlphaFold2...", flush=True)
        run_alphafold_local(msa_file, "output", local_db)
        return

    api_key = os.getenv('neurosnap_api_key')
    if not api_key:
        print("ERROR: API_KEY environment variable not set.", flush=True)
        sys.exit(1)

    # Extract the sequence text out of the generated fasta file
    sequence = ""
    with open(fasta_file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence += line.strip()

    print(f"Local DB not found. Preparing {id} for API endpoint submission...", flush=True)
    cleaned_sequence = "".join(sequence.split())
    input_sequence = {
        'aa': {id: cleaned_sequence},
        'rna': {},
        'dna': {}
    }

    submission_url = "https://neurosnap.ai/api/job/submit/AlphaFold2"

    multipart_data = MultipartEncoder(
        fields={
            "Input Sequences": json.dumps(input_sequence),
            "Model Type": "auto",
            "Template Mode": "none",
            "MSA Mode": "single_sequence", # Let Neurosnap handle MSA on their cluster due to firstline (>query) issues with msa file
            "Number Recycles": "5"
        }
    )

    headers = {"X-API-KEY": api_key, "Content-Type": multipart_data.content_type}
    print("Submitting job to API AlphaFold2 endpoint...", flush=True)
        
    request = requests.post(submission_url, headers=headers, data=multipart_data)
        
    # Explicitly print backend server errors before breaking on .json() parsing
    if request.status_code != 200:
        print(f"CRITICAL API ERROR ({request.status_code}): {request.text}", flush=True)
        sys.exit(1)
            
    job_id = request.json()

    POLLING_INTERVAL = 60
    print(f"Job submitted. Job ID: {job_id}. Polling for status every {POLLING_INTERVAL} seconds...", flush=True)

    while True:
        status_url = f"https://neurosnap.ai/api/job/status/{job_id}"
        try:
            status_res = requests.get(status_url, headers={"X-API-KEY": api_key})
            job_status = status_res.json().strip().lower()
        except Exception as err:
            print(f"Status check failed temporarily: {err}", flush=True)
            time.sleep(POLLING_INTERVAL)
            continue

        if job_status == "completed":
            job_result_url = f"https://neurosnap.ai/api/job/file/{job_id}/out/rank_1.pdb"
            result_response = requests.get(job_result_url, headers={"X-API-KEY": api_key})
            
            # Fallback pathing check if file location names shift on newer server updates
            if result_response.status_code != 200:
                job_result_url = f"https://neurosnap.ai/api/job/file/{job_id}/rank_1.pdb"
                result_response = requests.get(job_result_url, headers={"X-API-KEY": api_key})

            with open(f"{id}_alphafold_result.pdb", "wb") as f:
                f.write(result_response.content)
            print("Job completed successfully. Result downloaded.", flush=True)
            break
    
        elif job_status in ["failed", "deleted", "cancelled"]:
            print(f"Job status returned termination code: {job_status}. Exiting.", flush=True)
            sys.exit(1)
        
        elif job_status in ["running", "pending"]:
            print(f"Job status: {job_status}. Waiting {POLLING_INTERVAL} seconds...", flush=True)
            time.sleep(POLLING_INTERVAL)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--id', required=True)
    parser.add_argument('--msa', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--local_db', required=False, default='/mnt/c/databases/af2_db')
    args = parser.parse_args()
    
    run_alphafold(args.msa, args.id, args.local_db, args.fasta)