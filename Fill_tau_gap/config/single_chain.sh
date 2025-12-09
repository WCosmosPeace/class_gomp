#!/bin/bash
#SBATCH -p th02
#SBATCH -J 1        # 任务名
#SBATCH -o log_%j.out       # 标准输出
#SBATCH -e log_%j.err       # 错误输出
#SBATCH -N 1                # 1 个节点
#SBATCH -n 1              
#SBATCH --mem=64G           # 内存

ulimit -l unlimited

cobaya-run cmb_qdw_lens_bao.yaml --r
