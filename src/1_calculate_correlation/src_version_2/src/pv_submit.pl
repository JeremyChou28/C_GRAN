#!/bin/bash

#SBATCH --job-name=MSV000099738         # 作业名
#SBATCH --partition=cpu        # cpu 队列
#SBATCH -n 40 # 总核数 40
#SBATCH --ntasks-per-node=40    # 每节点核数
#SBATCH --exclusive

python -u entry.py \
  --intensity_file ../input/coral_processed_area/AREA_MSV000099738.txt \
  --compounds_num 82037 \
  --samples_num 375 \
  --correlation_result_filename ../output/coral_new/MSV000099738/correlation_results_MSV000099738_pos.csv \
  --n_jobs 8 \
  --skip_before_bh_csv \
  --tmp_name coral/MSV000099738_pos \
  > ../logs/coral_new/run_MSV000099738_pos_$(date +%F_%H-%M-%S).log 2>&1 