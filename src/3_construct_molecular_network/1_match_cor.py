import pandas as pd
from tqdm import tqdm

# 读取两个文件
file_a = pd.read_csv('input/source_target.csv')
file_b = pd.read_csv('output/corr_pval_final_CD_sediment_pos_3SD_20240828_miniCor_p0.05.csv')

# 转换为整数类型，确保匹配一致性
file_a['source'] = file_a['source'].astype(int)
file_a['target'] = file_a['target'].astype(int)
file_b['Substance 1'] = file_b['Substance 1'].astype(int)
file_b['Substance 2'] = file_b['Substance 2'].astype(int)

# 构建无序配对键（使 source/target 与 Substance 1/2 顺序无关）
file_a['pair_key'] = file_a[['source', 'target']].min(axis=1).astype(str) + '_' + file_a[['source', 'target']].max(axis=1).astype(str)
file_b['pair_key'] = file_b[['Substance 1', 'Substance 2']].min(axis=1).astype(str) + '_' + file_b[['Substance 1', 'Substance 2']].max(axis=1).astype(str)

# 精简 file_b 只保留用于合并的列
file_b_slim = file_b[['pair_key', 'Correlation', 'P-Value']]

# 合并数据并加入进度条
with tqdm(total=1, desc="Merging correlation data") as pbar:
    merged = pd.merge(file_a, file_b_slim, on='pair_key', how='left')
    pbar.update(1)

# 删除构造用的临时列
merged.drop(columns=['pair_key'], inplace=True)

# 保存输出文件
merged.to_csv('output/source_target_cor.csv', index=False)
