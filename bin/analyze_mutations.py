#!/usr/bin/env python3
import os
import sys
import json
import argparse
import numpy as np
from Bio import PDB
from Bio import SeqIO

def find_mutation_index_from_fasta(fasta_path, wt_header, mt_header):
    """Parses the FASTA file to find the 1-based index where WT and MT differ."""
    wt_seq = None
    mt_seq = None
    
    # Read FASTA entries dynamically using Biopython
    for record in SeqIO.parse(fasta_path, "fasta"):
        if record.id == wt_header:
            wt_seq = str(record.seq)
        elif record.id == mt_header:
            mt_seq = str(record.seq)
            
    if not wt_seq or not mt_seq:
        print(f"CRITICAL: Could not find sequences for {wt_header} or {mt_header} in FASTA.", file=sys.stderr)
        sys.exit(1)
        
    # Find the first mismatch
    min_len = min(len(wt_seq), len(mt_seq))
    for i in range(min_len):
        if wt_seq[i] != mt_seq[i]:
            return i + 1 # 1-based indexing for structural residues
            
    return None

def get_ca_atoms(structure):
    """Extracts all Alpha Carbon atoms from a structure."""
    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    ca_atoms.append(residue['CA'])
    return ca_atoms

def calculate_local_rmsd(wt_atoms, mt_atoms, center_res_idx, radius=8.0):
    """Calculates RMSD only for atoms within a spatial radius of a specific residue index."""
    target_ca = None
    for atom in wt_atoms:
        if atom.get_parent().get_id()[1] == center_res_idx:
            target_ca = atom
            break
            
    if target_ca is None:
        return 0.0

    center_coord = target_ca.get_coord()
    local_wt_atoms = []
    local_mt_atoms = []
    
    for wt_atom, mt_atom in zip(wt_atoms, mt_atoms):
        wt_coord = wt_atom.get_coord()
        distance = np.linalg.norm(wt_coord - center_coord)
        if distance <= radius:
            local_wt_atoms.append(wt_atom)
            local_mt_atoms.append(mt_atom)
            
    if not local_wt_atoms:
        return 0.0
        
    sq_distances = []
    for wt_a, mt_a in zip(local_wt_atoms, local_mt_atoms):
        diff = wt_a.get_coord() - mt_a.get_coord()
        sq_distances.append(np.dot(diff, diff))
        
    return np.sqrt(np.mean(sq_distances))

def extract_plddt(structure, res_idx):
    """Extracts pLDDT (stored in the B-factor column) for a specific residue."""
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[1] == res_idx:
                    if 'CA' in residue:
                        return residue['CA'].get_bfactor()
    return None

def main():
    parser = argparse.ArgumentParser(description="Structural Mutation Metrics Engine")
    parser.add_argument("--wt_pdb", required=True)
    parser.add_argument("--mt_pdb", required=True)
    parser.add_argument("--fasta", required=True, help="Path to the mutations template FASTA file")
    parser.add_argument("--wt_id", required=True, help="Header ID of WT sequence in FASTA (e.g., 2LZM_T4L_WT)")
    parser.add_argument("--mt_id", required=True, help="Header ID of Mutant sequence in FASTA (e.g., 2LZM_T4L_M99A)")
    parser.add_argument("--active_sites", type=str, default="", help="Space-separated functional residues")
    args = parser.parse_args()

    # 1. Dynamically compute mutation index inside Python
    mutation_idx = find_mutation_index_from_fasta(args.fasta, args.wt_id, args.mt_id)
    if mutation_idx == None:
        print(f"CRITICAL: No mutation mismatch found between {args.wt_id} and {args.mt_id}.", file=sys.stderr)
        sys.exit(1)

    # Parse space-separated active sites safely
    active_sites_list = [int(x) for x in args.active_sites.split()] if args.active_sites.strip() else []

    parser_pdb = PDB.PDBParser(QUIET=True)
    wt_struct = parser_pdb.get_structure("WT", args.wt_pdb)
    mt_struct = parser_pdb.get_structure("MT", args.mt_pdb)

    wt_ca = get_ca_atoms(wt_struct)
    mt_ca = get_ca_atoms(mt_struct)

    if len(wt_ca) != len(mt_ca):
        print("CRITICAL: Structural sequence lengths mismatch. Cannot align backbone reliably.", file=sys.stderr)
        sys.exit(1)

    # Global Alignment
    superimposer = PDB.Superimposer()
    superimposer.set_atoms(wt_ca, mt_ca)
    superimposer.apply(mt_struct.get_atoms())
    global_rmsd = superimposer.rms

    # Local Metrics calculation using our dynamic index
    mutation_local_rmsd = calculate_local_rmsd(wt_ca, mt_ca, mutation_idx, radius=8.0)

    active_site_rmsds = {}
    for site in active_sites_list:
        active_site_rmsds[f"site_{site}_rmsd"] = calculate_local_rmsd(wt_ca, mt_ca, site, radius=6.0)

    wt_plddt = extract_plddt(wt_struct, mutation_idx)
    mt_plddt = extract_plddt(mt_struct, mutation_idx)
    delta_plddt = (wt_plddt - mt_plddt) if (wt_plddt and mt_plddt) else None

    # Decision Logic
    verdict = "Neutral Variant"
    if global_rmsd > 3.5:
        verdict = "Loss of Function (Total Structural Misfolding)"
    elif mutation_local_rmsd > 2.0 or any(val > 2.0 for val in active_site_rmsds.values()):
        verdict = "Altered/Loss of Function (Active Site Distortion)"
    elif delta_plddt and delta_plddt > 25.0:
        verdict = "Loss of Function (Mutation-Induced Local Destabilization)"

    # Compile Results Block safely casting NumPy float32 to native Python floats
    output_metrics = {
        "calculated_mutation_index": int(mutation_idx),
        "global_rmsd": round(float(global_rmsd), 4),
        "mutation_site_local_rmsd": round(float(mutation_local_rmsd), 4),
        "active_site_deviations": {k: round(float(v), 4) for k, v in active_site_rmsds.items()},
        "wild_type_site_plddt": round(float(wt_plddt), 2) if wt_plddt is not None else None,
        "mutant_type_site_plddt": round(float(mt_plddt), 2) if mt_plddt is not None else None,
        "delta_plddt": round(float(delta_plddt), 2) if delta_plddt is not None else None,
        "final_classification": verdict
    }

    print(json.dumps(output_metrics, indent=4))

if __name__ == "__main__":
    main()