import glob
import os 
import random

def csv_prot_lig(prot_lig_dir,out_file):
    '''
    write csv with protein name, path to pdb, path to sdf of ligand 
    into a csv used as input into diffdockL for inference
    prot_lig_dir (str): name or path to directory with directories that have pdb of receptor and sdf of ligand files
    out_file (str): output csv with protein name, protein pdb path and ligand sdf path for diffdockL
    '''
    for prot in glob.glob(f"{prot_lig_dir}/**/*.pdb", recursive=True):
        complex_name = os.path.basename(os.path.dirname(prot))
        protein_path = f"data/{prot}"
        ligand_file = glob.glob(os.path.join(os.path.dirname(prot), "*.sdf"))[0]
        ligand_description = f'data/{ligand_file}'

        # create out_file file if dne
        open(out_file, 'w').write('complex_name,protein_path,ligand_description,protein_sequence\n') if not os.path.exists(out_file) else None

        # read and append info to our csv 
        with open(out_file, 'r+') as csv: 
            csv.seek(0) # move pointer to beginning to check header
            first_line = csv.readline().strip()  
            if first_line == 'complex_name,protein_path,ligand_description,protein_sequence':
                csv.seek(0, 2) # move pointer to end of file
                csv.write(f'{complex_name},{protein_path},{ligand_description},\n')
            else: 
                # should ony run the first time 
                csv.write('complex_name,protein_path,ligand_description,protein_sequence')
                csv.write(f'{complex_name},{protein_path},{ligand_description},\n')

def create_random_subset(file, subset_size, out_file):
    '''
    creates a smaller csv for subset of proteins from csv 
    file (str): csv with protein name, protein pdb path and ligand sdf path for diffdockL
    subset_size (int): how many proteins in the subset?
    out_file (str): csv output with appropriate header and the paths of subset prot & ligs 
    '''
    # read big csv
    with open(file, 'r') as file:
        lines = file.readlines()

    # choose random subset of 200 proteins 
    header = lines[0]
    data_lines = lines[1:]
    if len(data_lines) >= subset_size:
        random_subset = random.sample(data_lines, subset_size)
    else:
        # if fewer than 200 lines, select all
        print(f"The file has less than {subset_size} lines.")
        random_subset = data_lines  

    # write the random subset to a new csv
    with open(out_file, 'w', newline='') as csv_subset: 
        csv_subset.write(header)
        for line in random_subset:
            csv_subset.write(line)

def main():
    # create big csv if it dne
    big_csv='/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/cov_pdb_eval.csv'
    if not os.path.exists(big_csv):
        csv_prot_lig(prot_lig_dir='CovPDB_complexes', out_file=big_csv)
    # create subset 
    create_random_subset(file=big_csv, subset_size=200, out_file='subset_covpdb_eval.csv')

if __name__ == "__main__":
    main()


