import pandas as pd
import numpy as np
from tqdm import tqdm

# 读取第一个CSV文件
df1 = pd.read_csv('tmp/corr_pval_with_n_CD_sediment_pos_3SD_20240811-unfiltered.csv')
# 读取第二个CSV文件
df2 = pd.read_csv('tmp/corr_pval_with_n_CD_sediment_pos_3SD_20240811-filtered.csv')

# 检查两个DataFrame是否具有相同的结构
if df1.shape != df2.shape:
    print("两个文件的结构不一致，无法比较。")
    exit()

# 创建一个列表来存储不同的相关系数及P值
different_correlations = []

# 获取总行数并计算进度条更新的次数
total_rows = len(df1)
num_updates = total_rows // 100

# 创建进度条
progress_bar = tqdm(total=total_rows, desc="Comparing correlations")

# 记录不同的相关系数及其位置
counter = 0
for index, row in df1.iterrows():
    # 检查相同位置的相关系数是否相同
    corr1 = row['Correlation']
    corr2 = df2.iloc[index]['Correlation']
    pval1 = row['P-Value']
    pval2 = df2.iloc[index]['P-Value']
    n1 = row['n']
    n2 =df2.iloc[index]['n']
    if pd.isnull(corr1) and pd.isnull(corr2):
        continue  # 如果两个相关系数都是缺失值，则跳过比较
    elif pd.isnull(corr1) or pd.isnull(corr2):
        different_correlations.append([row['Substance 1'], row['Substance 2'], corr1, pval1, n1, corr2, pval2, n2])
    elif round(corr1, 9) != round(corr2, 9):
        different_correlations.append([row['Substance 1'], row['Substance 2'], corr1, pval1, n1, corr2, pval2, n2])
    
    counter += 1
    if counter % 100 == 0:
        progress_bar.update(100)

# 关闭进度条
progress_bar.close()

# 将不同的相关系数及P值写入CSV文件
different_correlations_df = pd.DataFrame(different_correlations, columns=['Substance 1', 'Substance 2', 'Correlation 1', 'P-Value 1', 'n1', 'Correlation 2', 'P-Value 2', 'n2'])
different_correlations_df.to_csv('tmp/Sediment_pos_3SD_20240812_different_correlations_with_n.csv', index=False)

print("已将不同的相关系数P值及n值写入到文件中。")
