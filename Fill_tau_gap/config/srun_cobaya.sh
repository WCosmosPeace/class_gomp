#!/bin/bash
#SBATCH -p th02
#SBATCH -J srun        # 任务名
#SBATCH -o log_%j.out       # 标准输出
#SBATCH -e log_%j.err       # 错误输出
#SBATCH -N 1                # 1 个节点
#SBATCH -n 4                # 4 个核
#SBATCH --mem=128G           # 内存

ulimit -l unlimited

srun -n 4 cobaya-run srun_cmb_qdw_lens_bao.yaml
