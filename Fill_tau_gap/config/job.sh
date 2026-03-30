#!/bin/bash
#SBATCH -p th02
#SBATCH -J FTG
#SBATCH -o log_%j.out
#SBATCH -e log_%j.err
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --mem=200G

ulimit -l unlimited

srun -n 5 cobaya-run cmb_qdw_lens_bao.yaml --r