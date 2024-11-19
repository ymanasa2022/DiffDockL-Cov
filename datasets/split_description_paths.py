protein_paths_file = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/PDBBind_processed/ligand_descriptions.txt'
subset_ids_file = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/splits/timesplit_pdbbind_subset_train'
output_file = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/PDBBind_processed/ligand_descriptions_train.txt'

# read list of protein paths
with open(protein_paths_file, 'r') as f:
    protein_paths = f.read().splitlines()

# extract protein IDs from the paths
protein_id_to_path = {path.split('/')[-2]: path for path in protein_paths}

# read the list of protein IDs 
with open(subset_ids_file, 'r') as f:
    subset_ids = f.read().splitlines()

# find the matching paths
matched_paths = [protein_id_to_path[protein_id] for protein_id in subset_ids if protein_id in protein_id_to_path]

# write the matched paths to an output file
with open(output_file, 'w') as f:
    for path in matched_paths:
        f.write(path + '\n')

print(f'Matched protein paths written to {output_file}')


