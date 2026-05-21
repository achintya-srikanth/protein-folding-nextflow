#!/usr/bin/env python3

import sys
import argparse
import requests
import time
import os
import subprocess
from requests_toolbelt.multipart.encoder import MultipartEncoder
import json

# Define parameters
parser = argparse.ArgumentParser()
parser.add_argument('--id', required=True)
parser.add_argument('--sequence', required=True)
parser.add_argument('--local_db', default="/mnt/c/databases/mmseqs_db")
args = parser.parse_args()

def run_mmseqs_local(fasta_path, db_path, output_a3m):
    # Standard ColabFold/MMseqs2 search workflow
    # 1. Search, 2. Align, 3. Generate A3M
    cmd = [
        "mmseqs colab5msearch",
        fasta_path,
        db_path,
        "tmp_dir",
        output_a3m
    ]
    subprocess.run(" ".join(cmd), shell=True, check=True)

def get_msa(sequence, id, local_db):
    fasta_path = f"{id}.fasta"
    with open(fasta_path, "w") as f:
        f.write(f">{id}\n{sequence}\n")

    if local_db and os.path.exists(local_db):
        print("Fetching MSA. Using local MMseqs2 database at {}".format(local_db))
        run_mmseqs_local(fasta_path, local_db, f"{id}_sequence_msa.a3m")
        return

    print("Fetching MSA. No local database found, using API Endpoint.")
    api_key = os.getenv('neurosnap_api_key')
    if not api_key:
        print("ERROR: NEUROSNAP_API_KEY environment variable is not set.")
        sys.exit(1)

    # 3. CONSTRUCT MULTIPART DATA
    url = "https://neurosnap.ai/api/job/submit/mmseqs2 MSA Generation"
    
    clean_sequence = "".join(sequence.split()) # Removes ALL whitespace, newlines
    input_sequences = {
        "aa": {id: clean_sequence},
        "dna": {},
        "rna": {}
    }

    print(f"Preparing MMseqs2 MSA Generation submission for: {id}")
    
    multipart_data = MultipartEncoder(
        fields={
            "Query Sequence": json.dumps(input_sequences),
            "Coverage": "80",          # Higher coverage for better alignments
            "Identity threshold": "30", # Standard for remote homologs
            "Max Sequences": "1000",
        }
    )

    headers = {
        "X-API-KEY": api_key,
        "Content-Type": multipart_data.content_type
    }

    # 4. SUBMIT JOB
    try:
        response = requests.post(url, headers=headers, data=multipart_data)
        response.raise_for_status()
        
        # Neurosnap returns the job ID as a raw JSON string
        job_id = response.json()
        print(f"Job submitted successfully. Job ID: {job_id}")
        
    except Exception as e:
        print(f"CRITICAL ERROR during submission: {e}")
        if response is not None:
            print(f"Server responded with: {response.text}")
        sys.exit(1)

    # 5. POLLING FOR RESULTS

    POLL_INTERVAL = 60  # Seconds between status checks
    status_url = f"https://neurosnap.ai/api/job/status/{job_id}"
    print(f"Starting polling for results (Interval: {POLL_INTERVAL}s)...")

    while True:
        try:
            status_res = requests.get(status_url, headers={"X-API-KEY": api_key})
            status_res.raise_for_status()
            
            # Neurosnap returns a raw string, not a dict!
            raw_status = status_res.json()
            job_status = str(raw_status).strip().lower() 
            
            print(f"DEBUG: Server says status is: {job_status}", flush=True)

        except Exception as poll_err:
            print(f"Status check error: {poll_err}", flush=True)
            time.sleep(POLL_INTERVAL)
            continue 

        if job_status in ["completed", "success", "done"]:
            print("Folding complete! Attempting download...", flush=True)
            
            # Using the exact URL structure you verified
            download_url = f"https://neurosnap.ai/api/job/file/{job_id}/out/final.a3m"
            
            try:
                pdb_res = requests.get(download_url, headers={"X-API-KEY": api_key})
                pdb_res.raise_for_status() 
                
                output_filename = f"{id}_sequence_msa.a3m"
                with open(output_filename, "wb") as f:
                    f.write(pdb_res.content)
                
                print(f"Success! MSA saved to {output_filename}", flush=True)
                break 
                
            except Exception as dl_err:
                print(f"DOWNLOAD FAILED: {dl_err}", flush=True)
                sys.exit(1) 
            
        elif job_status in ["failed", "error"]:
            print("ERROR: The Neurosnap job failed internally.", flush=True)
            sys.exit(1)
            
        else: 
            # Status is likely "running", "queued", or "submitted"
            time.sleep(POLL_INTERVAL)

if __name__ == "__main__":

    get_msa(args.sequence, args.id, args.local_db)
