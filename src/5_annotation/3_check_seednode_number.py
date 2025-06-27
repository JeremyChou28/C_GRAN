import pandas as pd
from tqdm import tqdm

# 输入输出路径
input_path = 'tmp/Seednode_and_Targetnode.csv'  # 替换为你的实际文件路径
output_unique = 'tmp/Seednode_and_Targetnode_unique_seednode.csv'
output_repeated = 'tmp/Seednode_and_Targetnode_not_unique_seednodes.csv'

# 读取CSV
df = pd.read_csv(input_path)

# 使用 tqdm 展示拆分进度
with tqdm(total=1, desc="Splitting based on Targetnode count") as pbar:
    # 计算每个 Targetnode 出现次数
    target_counts = df['Targetnode'].value_counts()

    # 获取只出现一次和多次的 Targetnode ID 集合
    unique_targets = set(target_counts[target_counts == 1].index)
    repeated_targets = set(target_counts[target_counts > 1].index)

    # 拆分为两个 DataFrame
    df_unique = df[df['Targetnode'].isin(unique_targets)].copy()
    df_repeated = df[df['Targetnode'].isin(repeated_targets)].copy()

    pbar.update(1)

# 保存两个文件
df_unique.to_csv(output_unique, index=False)
df_repeated.to_csv(output_repeated, index=False)
