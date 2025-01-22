import os 
import shutil 

list_files = 'pose_bust_works.txt'
dir_path = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/CSKDE95_datamol_af2'
output_dir = '/home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/data/working_pose_bust_cskde'

with open(list_files, 'r') as file:
    directories = [line.strip() for line in file if line.strip()] 

for dirs in directories: 
    src = os.path.join(dir_path, dirs)
    out = os.path.join(output_dir, dirs)

    shutil.copytree(src, out, dirs_exist_ok=True)
    