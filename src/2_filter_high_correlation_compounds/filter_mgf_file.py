import os
from tqdm import tqdm
import pandas as pd

# === 参数设置 ===
mgf_file_path = (
    "input/coral_spectrum_mgf_pos/"
    "GnpsMgf_0_202671937_99738.mgf"
)

csv_file_path = (
    "tmp/coral_ppcps_seednode/"
    "MSV000099738/high_correlated_nodes.csv"
)

output_mgf_path = (
    "output/coral_spectrum_mgf_pos_ppcps_p0.05/"
    "GnpsMgf_0_202671937_99738_ppcps_p0.05.mgf"
)

# === 读取 CSV 文件中的 ID ===
id_df = pd.read_csv(
    csv_file_path,
    dtype={"ID": "string"},
)

id_set = set(
    id_df["ID"]
    .dropna()
    .str.strip()
)

print(f"CSV 中读取到 {len(id_set)} 个唯一 ID")
print("ID 示例：", list(id_set)[:5])

# === 逐段读取并筛选 MGF ===
with open(mgf_file_path, "r") as f:
    lines = f.readlines()

selected_blocks = []
block = []
keep = False

# 统计保留的谱图数量
selected_spectrum_count = 0

for line in tqdm(lines, desc="筛选光谱块"):
    stripped = line.strip()

    if stripped == "BEGIN IONS":
        block = [line]
        keep = False

    elif stripped.startswith("SCANS="):
        scan_id = stripped.split("=", 1)[1].strip()

        if scan_id in id_set:
            keep = True

        block.append(line)

    elif stripped == "END IONS":
        block.append(line)
        block.append("\n")

        if keep:
            selected_blocks.extend(block)
            selected_spectrum_count += 1

        block = []

    else:
        block.append(line)

# === 创建输出目录 ===
output_dir = os.path.dirname(output_mgf_path)
if output_dir:
    os.makedirs(output_dir, exist_ok=True)

# === 写入输出文件 ===
with open(output_mgf_path, "w") as f:
    f.writelines(selected_blocks)

print(f"筛选完成，共保留 {len(selected_blocks)} 行")
print(f"最终 MGF 文件中共有 {selected_spectrum_count} 张谱图")
print(f"结果保存至：{output_mgf_path}")