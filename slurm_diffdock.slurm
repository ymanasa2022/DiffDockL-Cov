#!/bin/sh

#SBATCH --job-name=diffdock
#SBATCH --account=tromeara0
#SBATCH --partition=spgpu,gpu
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=20GB
#SBATCH --time=24:00:00
#SBATCH --output=diffdock_ex.out
#SBATCH --error=diffdock_ex.err

conda init
source ~/.bashrc
conda activate diffdock_L 
python -m inference --config default_inference_args.yaml  --protein_ligand_csv data/subset_covpdb_eval.csv --out_dir results/covpdb_200_inference
