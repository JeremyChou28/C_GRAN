import pandas as pd
from tqdm import tqdm

# 读取第一个CSV文件
file1 = pd.read_csv('corr_pval_with_n_CD_sediment_pos_3SD_20240811-unfiltered.csv')

# 读取第二个CSV文件
file2 = pd.read_csv('significant_Sediment_pos_3SD_20240828_different_correlations_with_n_true.csv')

# 使用Substance 1和Substance 2作为索引合并两个数据框
merged_df = pd.merge(file1, file2, on=['Substance 1', 'Substance 2'], how='left')

# 判断较小的Correlation值，并获取对应的P-Value
# 创建布尔掩码：True 表示 Correlation 1 更小，False 表示 Correlation 2 更小或相等
mask = merged_df['Correlation 1'] < merged_df['Correlation 2']

# 创建新的列 Correlation_min 和对应的 P-Value_min
merged_df['Correlation_min'] = merged_df['Correlation 1'].where(mask, merged_df['Correlation 2'])
merged_df['P-Value_min'] = merged_df['P-Value 1'].where(mask, merged_df['P-Value 2'])

# 将第一个文件中的 Correlation 和 P-Value 列替换为 Correlation_min 和 P-Value_min（若存在）
merged_df['Correlation'] = merged_df['Correlation_min'].fillna(merged_df['Correlation'])
merged_df['P-Value'] = merged_df['P-Value_min'].fillna(merged_df['P-Value'])

# 删除第二个文件中不需要的列
merged_df.drop(['Correlation 1', 'Correlation 2', 'P-Value 1', 'P-Value 2', 'Correlation_min', 'P-Value_min'], axis=1, inplace=True)

# 保存结果到新的CSV文件，只包含需要的四列，并添加进度条
with tqdm(total=merged_df.shape[0]) as pbar:
    merged_df[['Substance 1', 'Substance 2', 'Correlation', 'P-Value']].to_csv('corr_pval_final_CD_sediment_pos_3SD_20240828_miniCor.csv', index=False)
    pbar.update(merged_df.shape[0])
