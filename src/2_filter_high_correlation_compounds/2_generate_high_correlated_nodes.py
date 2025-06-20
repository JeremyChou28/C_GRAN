import pandas as pd
from tqdm import tqdm

# 读取原始CSV文件
df = pd.read_csv('corr_pval_final_CD_sediment_pos_3SD_20240828_miniCor_p0.05_connectwith_14seednode.csv')

# 筛选出 Correlation 大于 0.7 的行  (该阈值由用户自定义)
filtered_df = df[df['Correlation'] > 0.7]

# 初始化进度条：总共处理两列数据
with tqdm(total=2, desc="Processing IDs") as pbar:
    # 提取Substance 1和Substance 2列
    substances = pd.concat([filtered_df['Substance 1'], filtered_df['Substance 2']])
    pbar.update(1)

    # 去重并按升序排序
    unique_ids = pd.DataFrame(sorted(substances.unique()), columns=['ID'])
    pbar.update(1)

# 保存为新的CSV文件
unique_ids.to_csv('high_correlated_nodes.csv', index=False)
