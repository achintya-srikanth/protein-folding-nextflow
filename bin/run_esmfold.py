#!/usr/bin/env python3

import os
import requests
import argparse
import time
import sys
import subprocess
from requests_toolbelt.multipart.encoder import MultipartEncoder
import json

parser = argparse.ArgumentParser()
parser.add_argument('--id', required=True)
parser.add_argument('--sequence', required=True)
parser.add_argument('--local_db', required=False, default="/mnt/c/databases/esmfold_db")
args = parser.parse_args()

def run_esmfold_local(fasta_path, output_dir):
    # Documentation from EvolutionaryScale/esm
    cmd = [
        "esm-fold",
        f"-i {fasta_path}",
        f"-o {output_dir}/result.pdb"
    ]
    subprocess.run(" ".join(cmd), shell=True, check=True)

def run_esmfold(sequence, id, local_db):

    # Checking for local esmfold database for local execution
    fasta_path = f"{id}.fasta"
    with open(fasta_path, "w") as f:
        f.write(f">{id}\n{sequence}\n")

    if os.path.exists(local_db):
        print("Running ESMFold. Using local DB found at {}".format(local_db))
        output_dir = "esmfold_out"
        os.makedirs(output_dir, exist_ok=True)
        run_esmfold_local(fasta_path, output_dir)
        local_output = os.path.join(output_dir, "result.pdb")
        if os.path.exists(local_output):
            os.replace(local_output, f"{id}_esmfold_result.pdb")
        return

    # Local DB not found, using API endpoint
    print(f"Preparing ESMFold API Endpoint submission for: {id}")

    # Validate API Key
    api_key = os.getenv("neurosnap_api_key")
    if not api_key:
        print("ERROR: API Key not found. Please set the nextflow.config variable. Exiting...")
        sys.exit(1)
    
    print(f"Preparing sequence for {id} for submission...")

    cleaned_sequence = "".join(sequence.split()) # Remove whitespace/newlines
    input_sequence = {
        'aa': {id: cleaned_sequence},
        'rna': {},
        'dna': {}
    }


    multipart_data = MultipartEncoder(
        fields = {
            "Input Sequences": json.dumps(input_sequence),
            "Number Recycles": "6"
        }
    )

    headers = {
        "X-API-KEY": api_key,
        "Content-Type": multipart_data.content_type
    }

    # Make post request to submit job to ESMFold API with the sequence
    try:
        sequence_post_url = "https://neurosnap.ai/api/job/submit/ESMFold"
        request = requests.post(sequence_post_url, data=multipart_data, headers=headers)
        request.raise_for_status()
        job_id = request.json()

    except Exception as e:
        print(f"Exception Raised: {e}")
        if request is not None:
            print(f"Response Status Code: {request.status_code}")
            print(f"Response Content: {request.content}")
        sys.exit(1)

    POLL_INTERVAL = 60  # Seconds between status checks
    print(f"Job submitted successfully. Job ID: {job_id}. Polling for status every {POLL_INTERVAL} seconds...")
    
    while True:

        # Make a get request for status
        job_status_url = f"https://neurosnap.ai/api/job/status/{job_id}"
        status = requests.get(job_status_url, headers={"X-API-KEY": api_key}).json().strip().lower()

        # Check status and perform appropriate actions
        if status == 'completed':
            job_result_url = f"https://neurosnap.ai/api/job/file/{job_id}/out/{id}.pdb"
            pdb_res = requests.get(job_result_url, headers={"X-API-KEY": api_key})
            with open(f"{id}_esmfold_result.pdb", "wb") as f:
                f.write(pdb_res.content)
            break
            
        elif status == 'failed':
            print('Job status: failed')
            sys.exit(1)
            
        elif status == 'pending':
            print('Job status: pending, wait 60 seconds')
            time.sleep(POLL_INTERVAL)
            continue

        elif status == 'running':
            print('Job status: running, wait 60 seconds')
            time.sleep(POLL_INTERVAL)
            continue
            
        elif status == 'deleted':
            print('Job status: deleted, exiting...')
            sys.exit(1)

        elif status == 'canceled':
            print('Job status: canceled, exiting...')
            sys.exit(1)
            
        else:
            print('Job status: unknown, exiting')
            sys.exit(1)

if __name__ == "__main__":

    # calling the function to run esmfold, get the pdb file and write results
    run_esmfold(args.sequence, args.id, args.local_db)