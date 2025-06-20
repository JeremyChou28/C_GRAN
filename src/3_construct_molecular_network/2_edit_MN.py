import pandas as pd
from tqdm import tqdm

# 输入输出路径
input_file = 'output/source_target_cor.csv'  # 替换为你的实际文件名
output_file = 'output/source_target_cor_edit.csv'

# 读取CSV文件
df = pd.read_csv(input_file)

# 初始化进度条（2 步）
with tqdm(total=2, desc="Filtering data") as pbar:
    
    # Step 1: 保留 Correlation > 0.7 的行  #这个相关性的阈值为之前筛选高相关物质时用户输入的阈值
    df_filtered = df[df['Correlation'] > 0.7].copy()
    pbar.update(1)
    
    # Step 2: 删除 abs(source_RT - target_RT) < 0.01 且 Correlation > 0.9 的行
    condition = (df_filtered['Correlation'] > 0.9) & (abs(df_filtered['source_RT'] - df_filtered['target_RT']) < 0.01)
    df_final = df_filtered[~condition].copy()
    pbar.update(1)

# 保存结果到新文件
df_final.to_csv(output_file, index=False)
