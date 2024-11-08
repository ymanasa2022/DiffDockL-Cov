#!/bin/sh

#SBATCH --job-name=diffdockL_train
#SBATCH --account=tromeara0
#SBATCH --partition=spgpu,gpu
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=20GB
#SBATCH --time=5:00:00
#SBATCH --output=diffdockL_train.out
#SBATCH --error=diffdockL_train.err


conda init
conda activate diffdock_L
python -m train \
--project DiffDockL_train \
--run_name pdbbind_partial_score_model \
--test_sigma_intervals \
--pdbbind_esm_embeddings_path /home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/datasets/pdbbind_embeddings.pt \
--log_dir /home/ymanasa/turbo/ymanasa/opt/DiffDockL-Cov/ymanasa_train_test \
--cudnn_benchmark \
--n_epochs 10
