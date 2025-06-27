import os
import pandas as pd
from tqdm import tqdm

# 路径设置（请修改为你实际的路径）
file_a_path = 'tmp/Seednode_and_Targetnode_not_unique_seednodes.csv'  # 文件A路径
folder_b = 'tmp/Seednode_and_Targetnode_Morgan_Similarity_score_split_not_unique' # 文件夹B路径
output_folder = 'tmp/Seednode_and_Targetnode_Morgan_Similarity_score_split_not_unique_Top10' # 输出文件夹

os.makedirs(output_folder, exist_ok=True)

# 读取文件A，建立 (Targetnode, Seednode) → Correlation 映射
df_corr = pd.read_csv(file_a_path)
# df_corr['Seednode'] = df_corr['Seednode'].astype(int)
# df_corr['Targetnode'] = df_corr['Targetnode'].astype(str)  # 保留一位小数

# 构建字典索引
correlation_dict = {
    (int(row['Targetnode']), int(row['Seednode'])): row['Correlation']
    for _, row in df_corr.iterrows()
}

# 构建：Targetnode → 多文件 DataFrame 聚合表
targetnode_groups = {}

# 遍历文件夹B
csv_files = [f for f in os.listdir(folder_b) if f.endswith('.csv') and '_' in f]

with tqdm(total=len(csv_files), desc="Processing files") as pbar:
    for filename in csv_files:
        try:
            # 解析文件名：Targetnode_Seednode.csv
            base_name = filename[:-4]  # 去掉 .csv
            target_str, seed_str = base_name.split('_')
            targetnode = eval(target_str)
            seednode = eval(seed_str)

            # 读取文件并添加 Correlation 加权列
            file_path = os.path.join(folder_b, filename)
            df = pd.read_csv(file_path)
            if 'score' not in df.columns:
                tqdm.write(f"Skipped {filename}: no 'score' column.")
                pbar.update(1)
                continue

            # 获取对应 Correlation 值
            correlation = correlation_dict.get((targetnode, seednode))
            if correlation is None:
                tqdm.write(f"Skipped {filename}: no Correlation match.")
                pbar.update(1)
                continue

            # 加权打分
            df['weighted_score'] = df['score'] * correlation
            df['Targetnode'] = targetnode
            df['Seednode'] = seednode

            # 收集该 Targetnode 的所有行
            if targetnode not in targetnode_groups:
                targetnode_groups[targetnode] = []
            targetnode_groups[targetnode].append(df)

        except Exception as e:
            tqdm.write(f"Error processing {filename}: {str(e)}")
        pbar.update(1)

# 对每个 Targetnode 聚合并筛选 top-N weighted_score
for targetnode, df_list in tqdm(targetnode_groups.items(), desc="Filtering and saving"):
    combined_df = pd.concat(df_list, ignore_index=True)

    # 按 weighted_score 降序排列
    combined_df.sort_values(by='weighted_score', ascending=False, inplace=True)
    top_scores = combined_df['weighted_score'].tolist()

    # 处理并列：找第10名及相同分数
    if len(top_scores) <= 10:
        df_top = combined_df
    else:
        threshold = top_scores[9]
        df_top = combined_df[combined_df['weighted_score'] >= threshold]

    # 输出文件名：Targetnode 整数位.csv
    output_name = f"{int(float(targetnode))}.csv"
    output_path = os.path.join(output_folder, output_name)
    df_top.to_csv(output_path, index=False)
