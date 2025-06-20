#!/bin/bash

#SBATCH --job-name=test        # 作业名
#SBATCH --partition=cpu        # cpu 队列
#SBATCH -N 1                   # 1个节点
#SBATCH --ntasks-per-node=40   # 每节点核数
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err


#python also_calculate_n_text.py CD_sediment_pos.txt \
#               --need-transpose \
#               --n-jobs 40 \
#               --run-name with_n_CD_sediment_pos_3SD_20240811-unfiltered

#python also_calculate_n_text.py CD_sediment_pos.txt \
#                             --n-jobs 40 \
#                             --need-transpose \
#                             --filter-expt \
#                             --run-name with_n_CD_sediment_pos_3SD_20240811-filtered