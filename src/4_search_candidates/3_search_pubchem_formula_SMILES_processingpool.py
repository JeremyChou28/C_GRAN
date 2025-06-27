import os
import re
import gzip
import pandas as pd
from tqdm import tqdm  # 进度条
from multiprocessing import Pool

# 现有 CSV 结果文件所在文件夹
csv_input_folder = "tmp/search_pubchem_mw_result"
# TTL 文件夹路径（存储 SMILES 和 Formula 信息）
smiles_folder_path = "/home/zhoujianping/hjwrw/Emerging_Pollutant_Molecular_Networking/pubchem_database/download_data/smiles"
formula_folder_path = "/home/zhoujianping/hjwrw/Emerging_Pollutant_Molecular_Networking/pubchem_database/download_data/formula"

def process_ttl_smiles(file_path):
    data = []
    with gzip.open(file_path, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('@') or not line.strip():
                continue
            parts = line.strip().split()
            cid = parts[0].split(':')[-1].split('_')[0]
            smiles = eval(parts[2])
            data.append((cid, smiles))
    return data

def process_ttl_formula(file_path):
    data=[]
    with gzip.open(file_path, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('@') or not line.strip():
                continue
            parts = line.strip().split()
            cid=parts[0].split(':')[-1].split('_')[0]
            formula=eval(parts[2])
            data.append((cid, formula))
    return data

# 用于并行处理多个文件
def process_files_in_parallel(file_paths, process_function, num_processes=4):
    with Pool(processes=num_processes) as pool:
        results = list(tqdm(pool.imap(process_function, file_paths), total=len(file_paths), desc="Processing files", unit="file"))
    return results


smiles_files = [f for f in os.listdir(smiles_folder_path) if f.endswith('.ttl.gz')]
smiles_file_paths = [os.path.join(smiles_folder_path, file_name) for file_name in smiles_files]
print("Processing SMILES files in parallel...")
smiles_results = process_files_in_parallel(smiles_file_paths, process_ttl_smiles)
smiles_database = {}
for result in smiles_results:
    smiles_database.update(result)  

formula_files = [f for f in os.listdir(formula_folder_path) if f.endswith('.ttl.gz')]
formula_file_paths = [os.path.join(formula_folder_path, file_name) for file_name in formula_files]
print("Processing Formula files in parallel...")
formula_results = process_files_in_parallel(formula_file_paths, process_ttl_formula)
formula_database = {}
for result in formula_results:
    formula_database.update(result)  


# **处理单行数据**
def process_row(row):
    cid = str(row["CID"])
    smiles = smiles_database.get(cid, "")
    formula = formula_database.get(cid, "")
    row["SMILES"] = smiles
    row["Formula"] = formula
    return row

# **处理所有 CSV 文件**
output_folder = "tmp/search_pubchem_mw_result_row_level"
os.makedirs(output_folder, exist_ok=True)
csv_files = [f for f in os.listdir(csv_input_folder) if f.endswith(".csv")]

# **读取所有 CSV 文件并合并为一个 DataFrame**
all_data = []
for csv_file in csv_files:
    input_file_path = os.path.join(csv_input_folder, csv_file)
    df = pd.read_csv(input_file_path)
    df["SourceFile"] = csv_file  # 添加源文件信息以便后续拆分
    all_data.append(df)
combined_df = pd.concat(all_data, ignore_index=True)

# **确保 CID 列存在**
if "CID" not in combined_df.columns:
    print("⚠️ 缺少 CID 列，无法处理")
    exit()

# **统一 CID 数据格式**
combined_df["CID"] = combined_df["CID"].astype(str)

print("正在匹配 SMILES 和 Formula 信息...")

# **并行处理每一行**
with Pool(processes=16) as pool:
    processed_rows = list(tqdm(pool.imap(process_row, combined_df.to_dict("records")), total=len(combined_df), desc="处理行", unit="row"))

# **将处理后的数据转换回 DataFrame**
processed_df = pd.DataFrame(processed_rows)

# **拆分并保存到对应的 CSV 文件中**
for csv_file in csv_files:
    output_file_path = os.path.join(output_folder, csv_file)
    file_df = processed_df[processed_df["SourceFile"] == csv_file].drop(columns=["SourceFile"])
    file_df.to_csv(output_file_path, index=False)

print(f"✅ 匹配完成，结果保存在 {output_folder} 文件夹中")