import pandas as pd
import os
from tqdm import tqdm

# === 参数设置 ===
file_A = "input/unique_ID_MW_Ocean_OPEs_cor0.00001_adjusted_pvalue_cor_decoy.csv"  # 包含 ID 和 MW 的文件
file_B = "output/filtered_merge_all_edit_filtered.csv"  # 包含 monoisotopicmass 的主文件
output_folder = "output/candidates_Ocean_OPEs_cor0.00001_adjusted_pvalue_cor_decoy_5ppm"
os.makedirs(output_folder, exist_ok=True)

# === 读取数据 ===
df_a = pd.read_csv(file_A)
df_b = pd.read_csv(file_B)

# 检查列名存在
assert 'ID' in df_a.columns and 'MW' in df_a.columns, "文件A缺失 ID 或 MW 列"
assert 'monoisotopicmass' in df_b.columns, "文件B缺失 monoisotopicmass 列"

# === 对每个 ID 执行筛选并保存 ===
for _, row in tqdm(df_a.iterrows(), total=len(df_a), desc="筛选中"):
    id_ = str(int(row['ID'])).strip()  # 将ID转换为整数再转为字符串
    mw = float(row['MW'])

    # 计算 ppm 偏差筛选
    filtered = df_b[abs(df_b['monoisotopicmass'] - mw) / mw * 1e6 < 5]  #记得修改质量偏差！！！

    if not filtered.empty:
        # 重命名 canonicalsmiles 列为 SMILES（如果存在）
        if 'canonicalsmiles' in filtered.columns:
            filtered = filtered.rename(columns={'canonicalsmiles': 'SMILES'})

        output_path = os.path.join(output_folder, f"{id_}.csv")
        filtered.to_csv(output_path, index=False)

print(f"完成！结果保存在文件夹：{output_folder}")
