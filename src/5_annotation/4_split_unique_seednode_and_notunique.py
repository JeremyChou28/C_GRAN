import os
import shutil
import pandas as pd
from tqdm import tqdm

# 文件路径配置（请修改为你自己的路径）
unique_csv = 'output/Seednode_and_Targetnode_unique_seednode.csv'
not_unique_csv = 'output/Seednode_and_Targetnode_not_unique_seednodes.csv'
source_folder = 'output/Seednode_and_Targetnode_Morgan_Similarity_score'  # 包含所有8579.0_27.0格式文件的文件夹
output_folder_unique = 'output/Seednode_and_Targetnode_Morgan_Similarity_score_split_unique'
output_folder_not_unique = 'output/Seednode_and_Targetnode_Morgan_Similarity_score_split_not_unique'

# 创建输出文件夹
os.makedirs(output_folder_unique, exist_ok=True)
os.makedirs(output_folder_not_unique, exist_ok=True)

# 读取两个csv文件
df_unique = pd.read_csv(unique_csv)
df_not_unique = pd.read_csv(not_unique_csv)

# 构建命名集合（文件名：'Targetnode_Seednode.csv'）
def build_filename_set(df):
    return set(f"{row['Targetnode']}_{row['Seednode']}.csv" for _, row in df.iterrows())

unique_files = build_filename_set(df_unique)
not_unique_files = build_filename_set(df_not_unique)

# 遍历原始文件夹中的所有csv文件，按匹配规则分类拷贝
all_files = [f for f in os.listdir(source_folder) if f.endswith('.csv')]

with tqdm(total=len(all_files), desc="Classifying CSV files") as pbar:
    for file in all_files:
        source_path = os.path.join(source_folder, file)

        if file in unique_files:
            shutil.copy(source_path, os.path.join(output_folder_unique, file))
        elif file in not_unique_files:
            shutil.copy(source_path, os.path.join(output_folder_not_unique, file))
        # 忽略未匹配文件
        pbar.update(1)
