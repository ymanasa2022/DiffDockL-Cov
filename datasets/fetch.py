import requests 

def get_uniprot_id(pdb_id, pretty=False):
    '''
    returns Uniprot ID for a given PDB ID 
    '''
    full_url = "%s/%s/%s?pretty=%s" % ("https://www.ebi.ac.uk/pdbe/api", "/mappings/uniprot", pdb_id, str(pretty).lower())

    json_results = requests.get( full_url ).json() 

    # pull out the UniProt id. for this PDB id:
    uniprot_id = json_results[pdb_id] 
    uniprot_id2 = uniprot_id["UniProt"]
    uniprot_id3 = list(uniprot_id2.keys()) # a list of the UniProt ids
    # handeling multiple uniprots 
    if len(uniprot_id3) == 1:
        uniprot_id4 = uniprot_id3[0]
        return uniprot_id4
    elif len(uniprot_id3) > 1:
        print(f"Multiple UniProt IDs found for {pdb_id}: {uniprot_id3}")
    else:
        print(f"No UniProt ID found for: {pdb_id}")

def fetch_af2(pdb_id, output_path_file):
    '''
    get AlphaFold2 structure using PDB ID
    output_path_file: path to file *.pdb
    '''
    uniprot_id = get_uniprot_id(pdb_id)

    af2_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"

    response = requests.get(af2_url)
    if response.status_code == 200:
        with open(f"{output_path_file}", "wb") as f:
            f.write(response.content)
        print(f"AlphaFold2 structure for {uniprot_id} downloaded successfully!")
    else:
        print(f"Failed to download AlphaFold2 structure for {uniprot_id}.")
