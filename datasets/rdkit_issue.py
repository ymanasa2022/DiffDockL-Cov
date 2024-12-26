import subprocess
import os
import glob
from rdkit import Chem
from rdkit.Chem import rdmolops
import datamol as dm 
import prody as pr
from fetch import fetch_af2
from create_csv_input import csv_prot_lig

mol2_file = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/CSKDE95/Structures/CSKDE95/1qdq/lig-native-pre.mol2'
mol2_obabel_fixed = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/CSKDE95/Structures/CSKDE95/1qdq/lig-native-obabel.mol2'
complex_name = '1qdq'
sdf_file = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/CSKDE95/Structures/CSKDE95/1qdq/1qdq_ligand.sdf'
obabel_format_fix = f"obabel {mol2_file} -O {mol2_obabel_fixed}"
subprocess.run(obabel_format_fix, shell=True, check=True)
mol = Chem.MolFromMol2File(mol2_obabel_fixed ,sanitize=False, removeHs=False)
if mol:
    print('datamol fixing')
    mol_fix = dm.fix_mol(mol)
    if mol_fix:
        print('performing sani')
        mol_sani = dm.sanitize_mol(mol_fix)
        print('done sani')
    else:
        print(f'sanitization failed for {complex_name}_ligand')
        # err_lig.append(f'{complex_name}_ligand')
    if mol_sani:
        print('performing std')
        try:
            mol_std = dm.standardize_mol(mol_sani)
            print('done std')
        except Exception as e:
            print(e) 
        
    else:
        # err_lig.append(f'{complex_name}_ligand')
        print(f'standardization failed for {complex_name}_ligand')

    if mol_std:
        # resolve aromaticity issues
        rdmolops.Kekulize(mol_std, clearAromaticFlags=True)
        Chem.MolToMolFile(mol_std, sdf_file)
    else:
       print('mol loading failed with datamol, trying rdkit kekulize')