import pandas as pd
import numpy as np
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests

input_file = "tmp/coral/MSV000086120_pos/corr_pval_final_miniCor.csv"
output_file = "tmp/coral/MSV000086120_pos/corr_pval_final_miniCor_BH.csv"

print("Loading CSV...")
df = pd.read_csv(input_file)

print(f"Total rows: {len(df)}")

# 找到有效 p-value
valid_mask = df["P-Value"].notna()
valid_pvals = df.loc[valid_mask, "P-Value"].values

print(f"Valid p-values: {len(valid_pvals)}")

# 初始化列
df["adjusted p-value"] = np.nan

print("Running BH correction...")

# tqdm 只是用于显示阶段进度
with tqdm(total=3, desc="BH Correction") as pbar:

    # Step 1: BH 校正
    adjusted_pvals = multipletests(
        valid_pvals,
        method="fdr_bh"
    )[1]

    pbar.update(1)

    # Step 2: 写回 dataframe
    df.loc[valid_mask, "adjusted p-value"] = adjusted_pvals

    pbar.update(1)

    # Step 3: 保存文件
    df.to_csv(output_file, index=False)

    pbar.update(1)

print(f"Saved to: {output_file}")