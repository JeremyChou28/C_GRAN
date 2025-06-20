import os
import csv
import re
import pandas as pd
from tqdm import tqdm  # 进度条库

# 设置 CSV 输入文件路径
csv_input_file = "output/unique_ids_with_mw.csv"

# 设置 TTL 文件夹路径（请修改为你的TTL文件所在目录）
ttl_folder_path = "input/pubchem_2025_0218/MonoIsotopicWeight"

# 读取 CSV 输入文件
input_data = pd.read_csv(csv_input_file)

# 解析 CID 和 目标分子量
target_data = {}
for _, row in input_data.iterrows():
    cid = str(row[0]).strip()  # 第一列是物质 ID
    mw = float(row[1])  # 第二列是分子量
    target_data[cid] = mw

# 正则表达式匹配 TTL 文件中的 CID 和 分子量
cid_pattern = re.compile(r"descriptor:CID(\d+)_Mono_Isotopic_Weight")
mw_pattern = re.compile(r"sio:SIO_000300\s+([\d.]+)")

# 解析 TTL 文件
database = {}
ttl_files = [f for f in os.listdir(ttl_folder_path) if f.endswith(".ttl")]

print("正在解析 TTL 文件...")
for ttl_file in tqdm(ttl_files, desc="解析进度", unit="file"):
    file_path = os.path.join(ttl_folder_path, ttl_file)
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            cid_match = cid_pattern.search(line)
            mw_match = mw_pattern.search(line)
            if cid_match and mw_match:
                cid = cid_match.group(1)
                mw = float(mw_match.group(1))
                database[cid] = mw

# 计算偏差并筛选符合 2ppm 误差范围的物质  #误差范围由用户自定义
output_folder = "output/search_pubchem_mw_result"
os.makedirs(output_folder, exist_ok=True)

print("正在匹配 2ppm 误差范围的物质...")
for target_cid, target_mw in tqdm(target_data.items(), desc="匹配进度", unit="物质"):
    matched_cids = []
    
    for db_cid, db_mw in database.items():
        ppm_error = abs(db_mw - target_mw) / target_mw * 10**6
        if ppm_error < 2:
            matched_cids.append((db_cid, db_mw))

    # 如果有匹配项，保存为 CSV
    if matched_cids:
        output_file = os.path.join(output_folder, f"{target_cid}.csv")
        with open(output_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["CID", "Mono_Isotopic_Weight"])  # CSV 头部
            writer.writerows(matched_cids)

print(f"匹配完成，结果保存在 {output_folder} 文件夹中")
