import pandas as pd
from tqdm import tqdm

# 输入输出文件路径（请替换为你的实际路径）
input_file = 'annotation_results.csv'
output_file = 'annotation_results_filtered.csv'

# 读取CSV文件
df = pd.read_csv(input_file)

# 初始化进度条
with tqdm(total=2, desc="Filtering and selecting columns") as pbar:
    # Step 1: 筛选 CFM-ID_score > 0.7  #用户自定义阈值
    df_filtered = df[df['CFM-ID_score'] > 0.7].copy()
    pbar.update(1)

    # Step 2: 只保留指定列
    columns_to_keep = ['ID', 'CID', 'Mono_Isotopic_Weight', 'SMILES', 'Formula', 'CFM-ID_score']
    df_final = df_filtered[columns_to_keep]
    pbar.update(1)

# 保存结果
df_final.to_csv(output_file, index=False)
