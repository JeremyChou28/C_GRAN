import os
import gzip
import time
import argparse
from tqdm import tqdm
import numpy as np
import pickle as pk
from multiprocessing import Pool

def parse_args():
    parser = argparse.ArgumentParser(description="Processing pubchem database.")
    parser.add_argument(
        "--pubchem_mw_folder",
        type=str,
        required=True,
        help="pubchem molecular weight folder containing TTL files",
    )
    parser.add_argument(
        "--pubchem_smiles_folder",
        type=str,
        required=True,
        help="pubchem smiles folder containing TTL files",
    )
    parser.add_argument(
        "--pubchem_formula_folder",
        type=str,
        required=True,
        help="pubchem formula folder containing TTL files",
    )
    return parser.parse_args()


def process_ttl_mw(file_path):
    data = []
    with gzip.open(file_path, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('@') or not line.strip():
                continue
            parts = line.strip().split()
            cid = parts[0].split(':')[-1].split('_')[0]
            mw = eval(parts[2])
            data.append((cid, mw))
    return data

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


if __name__ == '__main__':
    start_time = time.time()
    args= parse_args()
    # 获取全量的mw数据
    mw_data_path = args.pubchem_mw_folder
    mw_files = [f for f in os.listdir(mw_data_path) if f.endswith('.ttl.gz')]
    # all_mw_data = []
    # for file_name in tqdm(mw_files, desc="Processing MW files"):
    #     mw_data = process_ttl_mw(os.path.join(mw_data_path, file_name))
    #     all_mw_data.extend(mw_data)
    mw_file_paths = [os.path.join(mw_data_path, file_name) for file_name in mw_files]
    all_mw_data = []
    print("Processing MW files in parallel...")
    mw_results = process_files_in_parallel(mw_file_paths, process_ttl_mw)
    for result in mw_results:
        all_mw_data.extend(result)    
    with open('./all_mw_data.pk', 'wb') as f:
        pk.dump(all_mw_data, f)
    print("len of all_mw_data:", len(all_mw_data))
    
    # 获取全量的smiles数据
    smiles_data_path = args.pubchem_smiles_folder
    smiles_files = [f for f in os.listdir(smiles_data_path) if f.endswith('.ttl.gz')]
    smiles_file_paths = [os.path.join(smiles_data_path, file_name) for file_name in smiles_files]
    all_smiles_data = []
    # for file_name in tqdm(smiles_files, desc="Processing SMILES files"):
    #     smiles_data = process_ttl_smiles(os.path.join(smiles_data_path, file_name))
    #     all_smiles_data.extend(smiles_data)
    print("Processing SMILES files in parallel...")
    smiles_results = process_files_in_parallel(smiles_file_paths, process_ttl_smiles)
    for result in smiles_results:
        all_smiles_data.extend(result)
    
    with open('./all_smiles_data.pk', 'wb') as f:
        pk.dump(all_smiles_data, f)
    print("len of all_smiles_data:", len(all_smiles_data))
    
    # 获取全量的formula数据
    formula_data_path = args.pubchem_formula_folder
    formula_files = [f for f in os.listdir(formula_data_path) if f.endswith('.ttl.gz')]
    formula_file_paths = [os.path.join(formula_data_path, file_name) for file_name in formula_files]
    all_formula_data = []
    # for file_name in tqdm(formula_files, desc="Processing formula files"):
    #     formula_data = process_ttl_formula(os.path.join(formula_data_path, file_name))
    #     all_formula_data.extend(formula_data)
    print("Processing Formula files in parallel...")
    formula_results = process_files_in_parallel(formula_file_paths, process_ttl_formula)
    for result in formula_results:
        all_formula_data.extend(result)
    with open('./all_formula_data.pk', 'wb') as f:
        pk.dump(all_formula_data, f)
    print("len of all_formula_data:", len(all_formula_data))
    print("All data processed and saved.")
        
    # all_mw_data = pk.load(open('./all_mw_data.pk', 'rb'))
    # all_smiles_data = pk.load(open('./all_smiles_data.pk', 'rb'))
    # all_formula_data = pk.load(open('./all_formula_data.pk', 'rb'))
    
    # 合并数据
    # 将 smiles和formula 数据预处理为字典（hash map）
    smiles_dict = dict(all_smiles_data)
    formula_dict = dict(all_formula_data)

    # 构造最终结果
    merged_data = {}
    for cid, mw in tqdm(all_mw_data, desc="Merging mw, smiles and formula"):
        merged_data[cid] = {'mw': mw, 'smiles': smiles_dict.get(cid), "formula": formula_dict.get(cid)}

    # 保存合并后的数据
    with open('./pubchem_database.pk', 'wb') as f:
        pk.dump(merged_data, f)
    print("Pubchem database saved successfully.")
    print("Spend time: ", time.time() - start_time)
    
