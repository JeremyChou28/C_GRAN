import os
import csv
import re
import gzip
import pandas as pd
from tqdm import tqdm  # 进度条库
from multiprocessing import Pool
# 设置 CSV 输入文件路径
csv_input_file = "tmp/unique_ID_MW.csv"

# 设置 TTL 文件夹路径（请修改为你的TTL文件所在目录）
ttl_folder_path = "/home/zhoujianping/hjwrw/Emerging_Pollutant_Molecular_Networking/pubchem_database/download_data/mw"

# 读取 CSV 输入文件
input_data = pd.read_csv(csv_input_file)

# 解析 CID 和 目标分子量
target_data = {}
for _, row in input_data.iterrows():
    cid = str(row[0]).strip()  # 第一列是物质 ID
    mw = float(row[1])  # 第二列是分子量
    target_data[cid] = mw

def process_ttl_mw(file_path):
    data = {}
    with gzip.open(file_path, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('@') or not line.strip():
                continue
            parts = line.strip().split()
            cid = parts[0].split(':')[-1].split('_')[0]
            mw = eval(parts[2])
            data[cid] = mw
    return data

# 用于并行处理多个文件
def process_files_in_parallel(file_paths, process_function, num_processes=4):
    with Pool(processes=num_processes) as pool:
        results = list(tqdm(pool.imap(process_function, file_paths), total=len(file_paths), desc="Processing files", unit="file"))
    return results


mw_files = [f for f in os.listdir(ttl_folder_path) if f.endswith('.ttl.gz')]
mw_file_paths = [os.path.join(ttl_folder_path, file_name) for file_name in mw_files]
print("Processing MW files in parallel...")
mw_results = process_files_in_parallel(mw_file_paths, process_ttl_mw)
database = {}
for result in mw_results:
    database.update(result)  

# 计算偏差并筛选符合 2ppm 误差范围的物质  #误差范围由用户自定义
output_folder = "tmp/search_pubchem_mw_result"
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
