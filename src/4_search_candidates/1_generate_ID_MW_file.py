import pandas as pd
from tqdm import tqdm

# 输入输出路径
input_file = 'test_files/source_target_cor_edit.csv'   # 替换为实际文件路径
output_file = 'tmp/unique_ID_MW.csv'

# 读取数据
df = pd.read_csv(input_file)

# 初始化进度条
with tqdm(total=2, desc="Extracting unique IDs and MWs") as pbar:
    # 提取 source 和 target 对应的 ID 和 MW
    source_df = df[['source', 'source_MW']].rename(columns={'source': 'ID', 'source_MW': 'MW'})
    target_df = df[['target', 'target_MW']].rename(columns={'target': 'ID', 'target_MW': 'MW'})
    pbar.update(1)
    
    # 合并并去重，按 ID 排序
    combined = pd.concat([source_df, target_df], ignore_index=True)
    unique_df = combined.drop_duplicates(subset='ID').sort_values(by='ID').reset_index(drop=True)
    pbar.update(1)

# 保存结果
unique_df.to_csv(output_file, index=False)
