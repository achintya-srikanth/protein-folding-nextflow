#!/usr/bin/env python3
import os
import requests
import argparse
import json
import time
import sys
import subprocess
from requests_toolbelt.multipart.encoder import MultipartEncoder


def run_rosettafold_local(fasta_path, output_dir, db_path):
    # Based on the RoseTTAFold2/run_RF2.sh documentation
    # Requires the database path to be exported to the environment
    os.environ['RF2_DB'] = db_path
    cmd = [
        "/app/RoseTTAFold2/run_RF2.sh",
        fasta_path,
        output_dir
    ]
    subprocess.run(" ".join(cmd), shell=True, check=True)

def run_rosettafold(id, sequence, local_db):
    fasta_path = f"{id}.fasta"
    with open(fasta_path, "w") as f:
        f.write(f">{id}\n{sequence}\n")

    if local_db and os.path.exists(local_db):
        print(f"DEBUG: Local RosettaFold2 DB found at {local_db}. Using local binary.")
        output_dir = "rosettafold_out"
        os.makedirs(output_dir, exist_ok=True)
        run_rosettafold_local(fasta_path, output_dir, local_db)
        local_output = os.path.join(output_dir, "result.pdb")
        if os.path.exists(local_output):
            os.replace(local_output, f"{id}_rosettafold_result.pdb")
        return

    # API Key Validation
    api_key = os.getenv('neurosnap_api_key')
    if not api_key:
        print("ERROR: NEUROSNAP_API_KEY environment variable is not set.")
        sys.exit(1)

    # CONSTRUCT MULTIPART DATA
    
    url = "https://neurosnap.ai/api/job/submit/RoseTTAFold2"
    
    clean_sequence = "".join(sequence.split()) # Removes ALL whitespace, newlines
    input_sequences = {
        "aa": {id: clean_sequence},
        "dna": {},
        "rna": {}
    }

    print(f"Preparing RoseTTAFold2 submission for: {id}")
    
    multipart_data = MultipartEncoder(
        fields={
            "Input Sequences": json.dumps(input_sequences),
            "Symmetry": "Unspecified / Unknown",
            "Order": "1",
            "MSA Method": "mmseqs2", # no custom MSA field for RF2, can't use prior .a3m file
            "Number Recycles": "6",
        }
    )

    headers = {
        "X-API-KEY": api_key,
        "Content-Type": multipart_data.content_type
    }

    # SUBMIT JOB
    try:
        response = requests.post(url, headers=headers, data=multipart_data)
        response.raise_for_status()
        
        # Neurosnap returns the job ID as a raw JSON string
        job_id = response.json()
        print(f"Job submitted successfully. Job ID: {job_id}")
        
    except Exception as e:
        print(f"ERROR during submission: {e}")
        if response is not None:
            print(f"Server responded with: {response.text}")
        sys.exit(1)

    # POLLING FOR RESULTS
    
    POLL_INTERVAL = 60  # Seconds between status checks

    status_url = f"https://neurosnap.ai/api/job/status/{job_id}"
    print(f"Starting polling for results (Interval: {POLL_INTERVAL}s)...")

    while True:
        try:
            status_res = requests.get(status_url, headers={"X-API-KEY": api_key})
            status_res.raise_for_status()
            
            # Neurosnap returns a raw string, not a dict
            raw_status = status_res.json()
            job_status = str(raw_status).strip().lower() 
            
            print(f"DEBUG: Server says status is: {job_status}", flush=True)

        except Exception as poll_err:
            print(f"Status check error: {poll_err}", flush=True)
            time.sleep(POLL_INTERVAL)
            continue 

        # --- Handle the Status ---
        if job_status in ["completed", "success", "done"]:
            print("Folding complete! Attempting download...", flush=True)
            
            # Using the exact URL structure you verified
            download_url = f"https://neurosnap.ai/api/job/file/{job_id}/out/output.pdb"
            
            try:
                pdb_res = requests.get(download_url, headers={"X-API-KEY": api_key})
                pdb_res.raise_for_status() 
                
                output_filename = f"{id}_rosettafold_result.pdb"
                with open(output_filename, "wb") as f:
                    f.write(pdb_res.content)
                
                print(f"Success! PDB saved to {output_filename}", flush=True)
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

    parser = argparse.ArgumentParser(description='Run RoseTTAFold2')
    parser.add_argument('--id', required=True, help='Protein ID')
    parser.add_argument('--sequence', required=True, help='Amino Acid Sequence')
    parser.add_argument('--local_db', required=False, default="/mnt/c/databases/rf2")
    
    args = parser.parse_args()
    run_rosettafold(args.id, args.sequence, args.local_db)