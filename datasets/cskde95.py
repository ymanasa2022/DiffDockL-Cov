import subprocess
import os
import glob

def process_structure_with_commands(prot_lig_dir, out_dir):
    '''
    preprocess the protein pdb files in /cskde95/*/ac/system-pre.pdb s.t. there are no LIG atoms,
    convert the mol2 ligand files in  /cskde95/*/ligand-native-pre.mol2 to sdf, 
    and move the files to appropriate pdb folder in output_dir
    
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

        # ligand processing
        mol2_path = f'{prot_lig_dir}/{complex_name}/lig-native-pre.mol2'
        sdf_path = f'{out_dir}/{complex_name}/{complex_name}_ligand.sdf'
        obabel_command = f"obabel {mol2_path} -O {sdf_path}"
        subprocess.run(obabel_command, shell=True, check=True)
        # check lig
        if not os.path.exists(mol2_path):
            err += complex_name
    return err 


def main():
    prot_lig_dir = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/CSKDE95/Structures/CSKDE95'
    out_dir = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/evaluation_CSKDE95/processed_structures'

    err = process_structure_with_commands(prot_lig_dir, out_dir)
    print(err)

if __name__ == "__main__":
    main()

