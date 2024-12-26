import subprocess
import os
import glob
from rdkit import Chem
from rdkit.Chem import rdmolops
import datamol as dm 
import prody as pr
from fetch import fetch_af2
from create_csv_input import csv_prot_lig

def process_structures(prot_lig_dir, out_dir, suffix = '_processed'):
    '''
    input:
    prot_lig_dir (str)- path to pdb files eg. /cskde95/*/ac/system-pre.pdb
    out_dir (str)- path to processed protein pdb s.t. there are no LIG atoms 
    suffix (str)- any suffix after pdb id in pdb file names eg. for a123, pdb file is: a123_processed.pdb, suffix = '_processed'
    output:
    err- list of proteins not processed
    '''
    err = []
    for prot in glob.glob(f"{prot_lig_dir}/*/ac/system-pre.pdb", recursive=True):
        print(f"processing {prot}")
        complex_name = os.path.basename(os.path.dirname(os.path.dirname(prot))) # pdb id
        print(complex_name)
        new_prot_dir = f'{out_dir}/{complex_name}'
        os.makedirs(new_prot_dir, exist_ok=True)
        new_prot_file = f'{new_prot_dir}/{complex_name}{suffix}.pdb'

        if not os.path.exists(new_prot_file):   # dont redo if exists 
            pdb = pr.parsePDB(prot)
            sequence = pdb.ca.getSequence()  # Extracts the sequence from CA atoms
            # fetch af2 
            if 'X' in sequence:
                print("Invalid sequence detected (contains 'X'). Fetching AF2 structure")
                if not os.path.exists(new_prot_file):
                    fetch_af2(pdb_id=complex_name, output_path_file=new_prot_file)

            # not af2 protein processing 
            else:
                if not os.path.exists(new_prot_file): # dont process if processed protein exists 
                    awk_command = f"awk '$4 != \"LIG\"' {prot} > {new_prot_file}"
                    subprocess.run(awk_command, shell=True, check=True)

        # check if protein was processed/fetched af2 
        if not os.path.exists(new_prot_file):
            err.append(complex_name)
    return err 

def process_ligands(prot_lig_dir, out_dir, tool='datamol'):
    '''
    prot_lig_dir: convert the mol2 ligand files in /cskde95/*/ligand-native-pre.mol2 to sdf, 
    out_dir: path to output of processed ligand sdf files  
    tool: either 'obabel' or 'rdkit' or 'datamol'. Default: datamol. pure obabel and rdkit conversion might cause inference issues with DiffDock

    err: list of ligands not converted to sdf from mol2
    '''
    err_lig = []
    for prot in glob.glob(f"{prot_lig_dir}/*/ac/system-pre.pdb", recursive=True):
        complex_name = os.path.basename(os.path.dirname(os.path.dirname(prot))) # pdb id
        print(f"processing {complex_name}_ligand")
        
        # ligand processing
        mol2_file = f'{prot_lig_dir}/{complex_name}/lig-native-post.mol2'
        mol2_obabel_fixed = f'{prot_lig_dir}/{complex_name}/lig-native-obabel.mol2'
        sdf_file = f'{out_dir}/{complex_name}/{complex_name}_ligand.sdf'

        if tool == 'rdkit':
            print('running rdkit for ligand conversion')
            # mol2 format fixing using obabel
            obabel_command = f"obabel {mol2_file} -O {mol2_obabel_fixed}"
            subprocess.run(obabel_command, shell=True, check=True)
            # rdkit mol2 to sdf 
            mol = Chem.MolFromMol2File(mol2_obabel_fixed, sanitize=False)
            Chem.MolToMolFile(mol, sdf_file, kekulize=False)

        elif tool == 'datamol':
            print('running datamol for ligand conversion')
            obabel_format_fix = f"obabel {mol2_file} -O {mol2_obabel_fixed}"
            subprocess.run(obabel_format_fix, shell=True, check=True)
            mol = Chem.MolFromMol2File(mol2_obabel_fixed ,sanitize=False, removeHs=False)
            if mol:
                mol_fix = dm.fix_mol(mol)
                if mol_fix:
                    mol_sani = dm.sanitize_mol(mol_fix)
                else:
                    print(f'sanitization failed for {complex_name}_ligand')
                    # err_lig.append(f'{complex_name}_ligand')
                if mol_sani:
                    mol_std = dm.standardize_mol(mol_sani)
                else:
                    # err_lig.append(f'{complex_name}_ligand')
                    print(f'standardization failed for {complex_name}_ligand')
                if mol_std:
                    # resolve aromaticity issues
                    rdmolops.Kekulize(mol_std, clearAromaticFlags=True)
                    Chem.MolToMolFile(mol_std, sdf_file) 
                else:
                    print('mol loading failed with datamol')

        else: # tool == 'obabel' or None
            print('running obabel for ligand conversion')
            obabel_command = f"obabel {mol2_file} -O {sdf_file}"
            subprocess.run(obabel_command, shell=True, check=True)

        # check if ligand was processed 
        if not os.path.exists(sdf_file) and f'{complex_name}_ligand' not in err_lig:
            err_lig.append(f'{complex_name}_ligand')
    return err_lig

def main():
    prot_lig_dir = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/CSKDE95/Structures/CSKDE95/'
    out_dir = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/evaluation_CSKDE95_datamol_af2/'
    
    tool = 'datamol'
    # # preprocessing structures/ligands
    # print("processing proteins...")
    # err = process_structures(prot_lig_dir, out_dir, tool)
    # print(err)
    print("processing ligands next...")
    err_lig = process_ligands(prot_lig_dir, out_dir, tool)
    print(err_lig)
    print(len(err_lig))

    # # csv for inference
    # out_file='/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/evaluation_CSKDE95_datamol_af2/evaluation_CSKDE95_datamol_af2.csv'
    # if not os.path.exists(out_file): 
    #     csv_prot_lig(prot_lig_dir=out_dir, out_file=out_file)
    # else: 
    #     print('CSV already exists. remove before creating new csv')

if __name__ == "__main__":
    main()

