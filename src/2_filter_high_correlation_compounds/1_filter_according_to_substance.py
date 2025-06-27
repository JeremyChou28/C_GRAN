import pandas as pd
from tqdm import tqdm

# 读取第一个CSV文件
file1 = pd.read_csv('../1_calculate_correlation/corr_pval_final_CD_sediment_pos_3SD_20240828_true_p0.05.csv')

# 读取第二个CSV文件
file2 = pd.read_csv('14seednode.csv')

# 提取第二个CSV文件中的ID列
substances = file2['ID'].tolist()

# 筛选Substance 1或者Substance 2在ID中的行
filtered_df = file1[(file1['Substance 1'].isin(substances)) | (file1['Substance 2'].isin(substances))]

# 保存筛选后的结果到新的CSV文件，添加进度条
with tqdm(total=filtered_df.shape[0]) as pbar:
    filtered_df.to_csv('corr_pval_final_CD_sediment_pos_3SD_20240828_miniCor_p0.05_connectwith_14seednode.csv', index=False)
    pbar.update(filtered_df.shape[0])
