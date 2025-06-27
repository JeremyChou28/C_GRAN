import pandas as pd
import numpy as np
from scipy.stats import norm
from tqdm import tqdm

def fisher_z_transform(r, n):
    # 检查相关系数是否在 (-1, 1) 区间内，避免除以零和取对数出现负数的情况
    valid_r = np.clip(r, -0.999999999, 0.999999999)
    return 0.5 * (pd.Series(r).apply(lambda x: 0.5 * (np.log(1 + x) - np.log(1 - x))) / np.sqrt(1 / (n - 3)))

def z_test(z1, z2, n1, n2):
    se = np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    z_statistic = (z1 - z2) / se
    p_value = 2 * (1 - norm.cdf(abs(z_statistic)))
    return p_value

# 读取CSV文件
df = pd.read_csv("tmp/Sediment_pos_3SD_20240812_different_correlations_with_n.csv")

# 初始化进度条
total_iterations = len(df)
with tqdm(total=total_iterations, desc="Processing Rows") as pbar:
    # 进行Fisher's z变换
    df['z1'] = np.where(df['n1'] <= 3, np.nan, fisher_z_transform(df['Correlation 1'], df['n1']))
    df['z2'] = np.where(df['n2'] <= 3, np.nan, fisher_z_transform(df['Correlation 2'], df['n2']))
    pbar.update(len(df))  # 更新进度条

    # 进行Z检验
    df['p_value'] = np.where((df['n1'] <= 3) | (df['n2'] <= 3), np.nan, z_test(df['z1'], df['z2'], df['n1'], df['n2']))
    pbar.update(len(df))  # 更新进度条

# 选择显著差异的行
significant_rows = df[(df['p_value'] < 0.05) | df['p_value'].isna()]

# 输出到CSV文件
significant_rows.to_csv("tmp/significant_Sediment_pos_3SD_20240828_different_correlations_with_n_true.csv", index=False)
