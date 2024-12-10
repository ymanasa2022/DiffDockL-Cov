import subprocess
import os
import glob
from rdkit import Chem
from create_csv_input import csv_prot_lig

def process_structure_with_commands(prot_lig_dir, out_dir, tool=None):
    '''
    out_dir: preprocess the protein pdb files in /cskde95/*/ac/system-pre.pdb s.t. there are no LIG atoms,
    prot_lig_dir: convert the mol2 ligand files in /cskde95/*/ligand-native-pre.mol2 to sdf, 
    and move the files to appropriate pdb folder in output_dir
    tool: either 'obabel' or 'rdkit'. Default: obabel. obabel conversion might cause inference issues with DiffDock

    err: ligands not converted to sdf from mol2
    '''
    err = []
    for prot in glob.glob(f"{prot_lig_dir}/*/ac/system-pre.pdb", recursive=True):
        print(prot)
        complex_name = os.path.basename(os.path.dirname(os.path.dirname(prot))) # pdb id

        # protein processing
        new_prot_dir = f'{out_dir}/{complex_name}'
        os.makedirs(new_prot_dir, exist_ok=True)
        new_prot_file = f'{new_prot_dir}/{complex_name}_processed.pdb'
        
        awk_command = f"awk '$4 != \"LIG\"' {prot} > {new_prot_file}"
        subprocess.run(awk_command, shell=True, check=True)

        # check if protein was processed 
        if not os.path.exists(new_prot_file):
            err.append(complex_name)

        # ligand processing
        mol2_file = f'{prot_lig_dir}/{complex_name}/lig-native-pre.mol2'
        sdf_file = f'{out_dir}/{complex_name}/{complex_name}_ligand.sdf'

        if tool == 'rdkit':
            print('running rdkit for conversion')
            mol = Chem.MolFromMol2File(mol2_file, sanitize=False)
            try:
                Chem.SanitizeMol(mol)
                writer = Chem.SDWriter(sdf_file)
                writer.write(mol)
                writer.close()

            except Exception as e:
                print('Failed to load/convert molecule using RdKit:', e)

        else: # tool == 'obabel' or None
            print('running obabel for conversion')
            obabel_command = f"obabel {mol2_file} -O {sdf_file}"
            subprocess.run(obabel_command, shell=True, check=True)

        # check if ligand was processed 
        if not os.path.exists(sdf_file):
            err.append(f'{complex_name}_ligand')

    return err 


def main():
    prot_lig_dir = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/CSKDE95/Structures/CSKDE95'
    out_dir = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/evaluation_CSKDE95_rdkit/processed_structures'
    
    tool = 'rdkit'
    # preprocessing structures
    if not os.path.exists(out_dir):
        err = process_structure_with_commands(prot_lig_dir, out_dir, tool='rdkit')
        print(err)
    else:
        print(f'Processed structures and ligands already exist here: {out_dir}. Delete before running cskde95.py')

    # csv for inference
    out_file='/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/evaluation_CSKDE95_rdkit/evaluation_CSKDE95_rdkit.csv'
    if not os.path.exists(out_file): 
        csv_prot_lig(prot_lig_dir=out_dir, out_file=out_file)
    else: 
        print('CSV already exists. remove before creating new csv')

if __name__ == "__main__":
    main()

