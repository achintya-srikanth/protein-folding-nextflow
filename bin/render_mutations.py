#!/usr/bin/env python3
import os
import sys
import json
import argparse
from pymol import cmd

def generate_structural_images(wt_pdb, mt_pdb, json_report, protein_id):
    """Automates structural alignment, coloring, and ray-tracing via PyMOL wrapper."""
    
    # 1. Load the dynamic mutation index computed by your analysis script
    with open(json_report, 'r') as f:
        metrics = json.load(f)
    mutation_idx = metrics["calculated_mutation_index"]

    # 2. Reset the PyMOL canvas session
    cmd.reinitialize()

    # 3. Load both structural PDB coordinates into distinct workspace objects
    wt_obj = "WT_Structure"
    mt_obj = "MT_Structure"
    cmd.load(wt_pdb, wt_obj)
    cmd.load(mt_pdb, mt_obj)

    # 4. Perform global backbone structural alignment
    print(f"Aligning {mt_obj} to {wt_obj} globally...")
    alignment_summary = cmd.align(f"{mt_obj} & name CA", f"{wt_obj} & name CA")

    # 5. Color the scaffolds
    cmd.color("gray70", wt_obj)      # Neutral gray background for Wild Type
    cmd.color("aquamarine", mt_obj)  # Soft blue/green for Mutant tracing

    # 6. Isolate and focus on the mutation site microenvironment
    mutation_sphere_selection = "mutation_site_neighborhood"
    cmd.select(mutation_sphere_selection, f"byres ({wt_obj} & resi {mutation_idx} expand 8.0)")
    cmd.show("ribbon", mutation_sphere_selection)
    
    # 7. Render the specific mutated residue as explicit chemical sticks
    wt_stick_selection = f"{wt_obj} & resi {mutation_idx}"
    mt_stick_selection = f"{mt_obj} & resi {mutation_idx}"
    cmd.show("sticks", wt_stick_selection)
    cmd.show("sticks", mt_stick_selection)
    
    cmd.color("salmon", wt_stick_selection)
    cmd.color("yellow", mt_stick_selection)

    # 8. Configure camera view and render parameters
    cmd.orient(mutation_sphere_selection)
    cmd.zoom(mutation_sphere_selection, buffer=10.0)
    cmd.set("ray_opaque_background", 0)
    cmd.set("cartoon_fancy_helices", 1)

    # 9. Render image file to disk
    output_image_path = f"{protein_id}_mutation_overlay.png"
    cmd.png(output_image_path, width=1200, height=1200, dpi=300, ray=1)
    cmd.delete("all")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automated PyMOL Rendering Engine")
    parser.add_argument("--wt_pdb", required=True)
    parser.add_argument("--mt_pdb", required=True)
    parser.add_argument("--json_report", required=True, help="Path to json report containing calculated mutation index")
    parser.add_argument("--id", required=True)
    
    args = parser.parse_args()
    generate_structural_images(args.wt_pdb, args.mt_pdb, args.json_report, args.id)