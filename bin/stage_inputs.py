#!/usr/bin/env python3
import os
import json
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Dynamic PDB Pairing Engine")
    parser.add_argument("--pdb_dir", default="results", help="Root results directory")
    parser.add_argument("--out_json", default="staged_pairs.json", help="Output path for staged channel data")
    args = parser.parse_args()

    root_dir = Path(args.pdb_dir)
    all_files = list(root_dir.glob("*/*_result.pdb"))
    
    # Structure to hold found files organized by: [protein_group][model][condition] = file_path
    # Example: database['1BTL_TEM1']['alphafold']['E166A'] = Path(...)
    database = {}

    for p in all_files:
        filename = p.name
        # Remove extension
        base_name = filename.replace("_result.pdb", "")
        tokens = base_name.split("_")
        
        if len(tokens) < 4:
            continue # Skip files that don't match our naming convention
            
        pdb_id = tokens[0]
        prot_name = tokens[1]
        group_key = f"{pdb_id}_{prot_name}"
        
        model = tokens[-1]
        cond = "_".join(tokens[2:-1]) # Handles variants like E166A or complex mutant names safely
        
        if group_key not in database:
            database[group_key] = {}
        if model not in database[group_key]:
            database[group_key][model] = {}
            
        database[group_key][model][cond] = str(p.resolve())

    # Now, find valid pairs (Every mutant must have a corresponding WT for that model)
    staged_pairs = []
    
    for group_key, models in database.items():
        for model, conditions in models.items():
            if "WT" not in conditions:
                print(f"Warning: Found structural files for {group_key} ({model}) but missing the WT reference structure. Skipping.")
                continue
                
            wt_path = conditions["WT"]
            
            for cond, mt_path in conditions.items():
                if cond == "WT":
                    continue
                
                # Extract the base 4-letter PDB code for the UniProt API node map
                pdb_code = group_key.split("_")[0]
                
                pair_entry = {
                    "pdb_id": pdb_code,
                    "model": model,
                    "wt_path": wt_path,
                    "mt_path": mt_path,
                    "meta_id": f"{group_key}_{cond}_{model}",
                    "wt_id": f"{group_key}_WT",
                    "mt_id": f"{group_key}_{cond}"
                }
                staged_pairs.append(pair_entry)

    # Save the paired targets out to a clean JSON manifest
    with open(args.out_json, "w") as f:
        json.dump(staged_pairs, f, indent=4)
        
    print(f"Successfully staged {len(staged_pairs)} unique multi-model evaluation pairs inside {args.out_json}")

if __name__ == "__main__":
    main()