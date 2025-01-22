import os 
import shutil 
from Bio.PDB import PDBList 

list_files = 'pose_bust_works.out'
# dir_path = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/CSKDE95_datamol_af2'
output_dir = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/work_posebust_csdke_orig'

with open(list_files, 'r') as file:
    pdb_codes = [line.strip() for line in file if line.strip()] 


pdb_list = PDBList()

for pdb_code in pdb_codes:
    os.makedirs(os.path.join(output_dir, pdb_code), exist_ok=True)
    pdir=os.path.join(output_dir, pdb_code)
    pdb_list.retrieve_pdb_file(pdb_code=pdb_code, file_format='pdb', pdir=pdir)
    os.rename(os.path.join(pdir, f'pdb{pdb_code}.ent'), os.path.join(pdir, f'pdb{pdb_code}.pdb'))
