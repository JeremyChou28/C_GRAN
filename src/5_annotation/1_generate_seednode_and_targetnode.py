import pandas as pd 
from tqdm import tqdm

# 输入路径
file_a_path = 'input/source_target_cor_edit.csv'  # 替换为实际路径
file_b_path = 'input/14seednode.csv'              # 替换为实际路径
output_path = 'output/Seednode_and_Targetnode.csv'

# 读取文件
df_a = pd.read_csv(file_a_path)
df_b = pd.read_csv(file_b_path)

# 提取种子物质 ID 集合
seed_ids = set(df_b['ID'].astype(int))  # 确保为整数类型以匹配

# 初始化列表用于存储结果
results = []

# 使用 tqdm 包装 DataFrame 行迭代器
for _, row in tqdm(df_a.iterrows(), total=len(df_a), desc="Filtering connected nodes"):
    source = int(row['source'])
    target = int(row['target'])
    corr = row['Correlation']  # 提取 Correlation 值

    source_in_seed = source in seed_ids
    target_in_seed = target in seed_ids

    # 仅保留一端在种子集合中，另一端不在
    if source_in_seed and not target_in_seed:
        results.append({
            'Seednode': source,
            'Targetnode': f"{target:.1f}",
            'Correlation': corr
        })
    elif target_in_seed and not source_in_seed:
        results.append({
            'Seednode': target,
            'Targetnode': f"{source:.1f}",
            'Correlation': corr
        })

# 转换为 DataFrame 并按 Seednode 升序排序
result_df = pd.DataFrame(results)
result_df.sort_values(by='Seednode', inplace=True)

# 保存输出文件
result_df.to_csv(output_path, index=False)
