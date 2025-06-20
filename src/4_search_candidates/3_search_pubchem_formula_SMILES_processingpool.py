import os
import re
import pandas as pd
from tqdm import tqdm  # 进度条
from multiprocessing import Pool

# 现有 CSV 结果文件所在文件夹
csv_input_folder = "output/search_pubchem_mw_result"
# TTL 文件夹路径（存储 SMILES 和 Formula 信息）
smiles_folder_path = "input/pubchem_2025_0218/canSMILES"
formula_folder_path = "input/pubchem_2025_0218/formula"
# **修正正则表达式，严格匹配 Canonical_SMILES 和 Molecular_Formula**
cid_smiles_pattern = re.compile(r"descriptor:CID(\d+)_Canonical_SMILES")
cid_formula_pattern = re.compile(r"descriptor:CID(\d+)_Molecular_Formula")
data_pattern = re.compile(r'sio:SIO_000300\s+"([^"]+)"')  # 提取双引号中的数据

# **解析 TTL 文件，构建 CID -> SMILES / Formula 映射**
def parse_ttl_files(folder_path, data_type="SMILES"):
    database = {}
    ttl_files = [f for f in os.listdir(folder_path) if f.endswith(".ttl")]
    print(f"正在解析 {data_type} TTL 文件...")
    for ttl_file in tqdm(ttl_files, desc=f"解析 {data_type}", unit="file"):
        file_path = os.path.join(folder_path, ttl_file)
        with open(file_path, "r", encoding="utf-8") as f:
            for line in f:
                # 根据 data_type 选择不同的正则匹配
                cid_match = cid_smiles_pattern.search(line) if data_type == "SMILES" else cid_formula_pattern.search(line)
                data_match = data_pattern.search(line)
                if cid_match and data_match:
                    cid = cid_match.group(1)  # 提取 CID
                    value = data_match.group(1)  # 提取 SMILES 或 Formula
                    database[cid] = value
    return database

# **解析 SMILES 和 Formula 数据**
smiles_database = parse_ttl_files(smiles_folder_path, "SMILES")
formula_database = parse_ttl_files(formula_folder_path, "Formula")

# **处理单行数据**
def process_row(row):
    cid = str(row["CID"])
    smiles = smiles_database.get(cid, "")
    formula = formula_database.get(cid, "")
    row["SMILES"] = smiles
    row["Formula"] = formula
    return row

# **处理所有 CSV 文件**
output_folder = "output/search_pubchem_mw_result_row_level"
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