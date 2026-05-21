#!/usr/bin/env python3
import sys
import requests

def get_uniprot_id_from_pdb(pdb_id):
    """Queries UniProt to find the Accession ID associated with a PDB structure."""
    url = f"https://rest.uniprot.org/uniprotkb/search?query=database:pdb-{pdb_id}&fields=accession"
    response = requests.get(url, headers={"Accept": "application/json"})
    
    if response.status_code == 200:
        results = response.json().get('results', [])
        if results:
            return results[0].get('primaryAccession')
    return None

def fetch_active_sites(uniprot_id):
    """Fetches catalytic active sites using the Proteins API."""
    url = f"https://upward.ebi.ac.uk/proteins/api/features/{uniprot_id}?type=ACT_SITE"
    response = requests.get(url, headers={"Accept": "application/json"})
    
    if response.status_code != 200:
        return []
        
    features = response.json().get('features', [])
    positions = [int(feat['begin']) for feat in features]
    return positions

if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_id = sys.argv[1].upper()
        
        # If it looks like a 4-character PDB ID, map it to UniProt first
        if len(input_id) == 4:
            uniprot_id = get_uniprot_id_from_pdb(input_id)
        else:
            uniprot_id = input_id
            
        if uniprot_id:
            sites = fetch_active_sites(uniprot_id)
            print(" ".join(map(str, sites)))
        else:
            # Output nothing if mapping fails, downstream process handles empty string gracefully
            print("")