#!/bin/sh

#SBATCH --job-name=diffdock
#SBATCH --account=maom99
#SBATCH --partition=spgpu,gpu
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=20GB
#SBATCH --time=5:00:00
#SBATCH --output=diffdock_ex.out
#SBATCH --error=diffdock_ex.err


python -m inference --config /home/ymanasa/turbo/ymanasa/opt/DiffDock/default_inference_args.yaml --protein_ligand_csv data/protein_ligand_example.csv --out_dir results/user_predictions_small


