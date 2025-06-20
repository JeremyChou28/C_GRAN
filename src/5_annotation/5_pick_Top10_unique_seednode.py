import os
import pandas as pd
from tqdm import tqdm

# 输入与输出文件夹路径（请根据实际路径修改）
input_folder = 'output/Seednode_and_Targetnode_Morgan_Similarity_score_split_unique'       # 原始csv文件所在文件夹
output_folder = 'output/Seednode_and_Targetnode_Morgan_Similarity_score_split_unique_Top10' # 输出文件夹

# 创建输出文件夹
os.makedirs(output_folder, exist_ok=True)

# 获取所有以 "_xxx.csv" 命名的csv文件
csv_files = [f for f in os.listdir(input_folder) if f.endswith('.csv') and '_' in f]

# 初始化进度条
with tqdm(total=len(csv_files), desc="Filtering top score rows") as pbar:
    for filename in csv_files:
        filepath = os.path.join(input_folder, filename)
        try:
            df = pd.read_csv(filepath)
            if 'score' not in df.columns:
                tqdm.write(f"Skipped {filename}: no 'score' column")
                pbar.update(1)
                continue

            # 获取前10名及并列的 score
            scores = df['score'].tolist()
            if len(scores) <= 10:
                df_filtered = df
            else:
                threshold_score = sorted(scores, reverse=True)[9]
                df_filtered = df[df['score'] >= threshold_score]

            # 输出文件名：取第一个数值部分整数形式
            prefix = filename.split('_')[0]       # 例如 '262.0'
            output_filename = f"{int(float(prefix))}.csv"
            output_path = os.path.join(output_folder, output_filename)

            # 保存（覆盖）
            df_filtered.to_csv(output_path, index=False)

        except Exception as e:
            tqdm.write(f"Error processing {filename}: {str(e)}")
        pbar.update(1)
