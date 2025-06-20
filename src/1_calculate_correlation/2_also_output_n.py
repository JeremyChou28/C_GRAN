import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm

def load_pickle(pkl_path):
    f = open(pkl_path, 'rb')
    pkl = pickle.load(f)
    return pkl['corr'], pkl['p_value'], pkl['non_zero_count']

filtered_pkl_path = "corr_pval_with_n_CD_sediment_pos_3SD_20240811-filtered.pickle"
f_corr_df, f_pv_df , f_n_df = load_pickle(filtered_pkl_path)

# 获取物质名称
substances = f_corr_df.index.tolist()

# 创建一个空的DataFrame用于存储数据
combined_data = []

# 计算总共需要循环的次数
total_iterations = sum(range(len(substances)))

# 使用tqdm来显示进度条
with tqdm(total=total_iterations, desc="Progress") as pbar:
    # 循环遍历所有物质对，并将相关系数和P值和n值存储到新的DataFrame中
    for i, substance1 in enumerate(substances):
        for j, substance2 in enumerate(substances):
            if j < i:  # 只输出下三角部分
                correlation = f_corr_df.iloc[i, j]
                p_value = f_pv_df.iloc[i, j]
                n = f_n_df.iloc[i,j]
                combined_data.append([substance1, substance2, correlation, p_value, n])
                if (i * len(substances) + j) % 100 == 0:  # 每100次更新一次进度
                    pbar.update(100)  # 更新进度条

# 创建DataFrame并设置列名
combined_df = pd.DataFrame(combined_data, columns=['Substance 1', 'Substance 2', 'Correlation', 'P-Value', 'n'])

# 将数据保存到CSV文件中
combined_df.to_csv('corr_pval_with_n_CD_sediment_pos_3SD_20240811-filtered.csv', index=False)
