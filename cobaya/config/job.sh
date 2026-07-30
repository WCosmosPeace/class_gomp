#!/bin/bash
#SBATCH -J 1GeV
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

#SBATCH -N 1
#SBATCH --ntasks=15
#SBATCH --cpus-per-task=16

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

ulimit -l unlimited
mkdir -p logs

srun --cpu-bind=cores cobaya-run cmb_qdw_lens_bao.yaml --resume
