import os
import csv
import re
import gzip
import time
import argparse
import pandas as pd
import pickle as pk
import pandas as pd
import bisect
from tqdm import tqdm
from multiprocessing import Pool


def parse_args():
    parser = argparse.ArgumentParser(description="Search candidates.")
    # online parameters
    parser.add_argument(
        "--ppm_threshold",
        default=2,
        type=int,
        required=True,
        help="the ppm threshold for searching candidates based on molecular weight",
    )
    parser.add_argument(
        "--is_filter_element",
        action="store_true",
        default=True,
        help="whether to filter candidates by element set",
    )
    parser.add_argument(
        "--element_set",
        type=str,
        required=True,
        help="the element set to filter candidates, e.g., 'C,H,O,N,P,S,F,Cl,Br,I'",
    )
    # offline parameters
    parser.add_argument(
        "--edited_molecular_network_file",
        default="source_target_cor_edit.csv",
        type=str,
        required=True,
        help="input molecular network file containing source, target, corrrelation",
    )
    parser.add_argument(
        "--pubchem_database_path",
        type=str,
        required=True,
        help="pubchem database path",
    )
    parser.add_argument(
        "--candidates_folder",
        type=str,
        required=True,
        help="the candidates folder to save the search results",
    )
    return parser.parse_args()


def generate_ID_MW_file(molecular_network_df, output_filename):
    # 初始化进度条
    with tqdm(total=2, desc="Extracting unique IDs and MWs") as pbar:
        # 提取 source 和 target 对应的 ID 和 MW
        source_df = molecular_network_df[["source", "source_MW"]].rename(
            columns={"source": "ID", "source_MW": "MW"}
        )
        target_df = molecular_network_df[["target", "target_MW"]].rename(
            columns={"target": "ID", "target_MW": "MW"}
        )
        pbar.update(1)

        # 合并并去重，按 ID 排序
        combined = pd.concat([source_df, target_df], ignore_index=True)
        unique_df = (
            combined.drop_duplicates(subset="ID")
            .sort_values(by="ID")
            .reset_index(drop=True)
        )
        pbar.update(1)

    # 保存结果
    unique_df.to_csv(output_filename, index=False)
    return unique_df


def find_candidates_by_ppm_smiles(target_mws, pubchem_mws, ppm_threshold):
    """
    target_mws: list of (target_id, mw)
    pubchem_mws: list of (cid, mw, smiles, formula)
    return: dict, key is cid, value is list of (cid, mw, smiles, formula)
    """
    print("Sorting pubchem_mws by molecular weight...")
    pubchem_mws = sorted(pubchem_mws, key=lambda x: x[1])

    pubchem_mw_values = [mw for _, mw, _, _ in tqdm(pubchem_mws, desc="Extracting MWs")]
    # pubchem_cid_smiles_formula = [(cid, smiles, formula) for cid, _, smiles, formula in tqdm(pubchem_mws, desc="Extracting CID & SMILES & Formula")]

    result = {}
    for target_id, target_mw in tqdm(target_mws, desc="Finding candidate SMILES"):
        delta = target_mw * ppm_threshold * 1e-6
        lower = target_mw - delta
        upper = target_mw + delta

        left = bisect.bisect_left(pubchem_mw_values, lower)
        right = bisect.bisect_right(pubchem_mw_values, upper)

        candidates = pubchem_mws[left:right]

        if candidates:
            result[target_id] = candidates

    return result


def save_results_to_csv(result_dict, save_dir="./result_csv"):
    """
    result_dict: dict[target_id] = list of (cid, smiles, formula)
    Save each target_id's result into one CSV file.
    """
    os.makedirs(save_dir, exist_ok=True)

    for target_id, candidates in tqdm(result_dict.items(), desc="Saving CSV files"):
        save_path = os.path.join(save_dir, f"{target_id}.csv")
        with open(save_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["CID", "MW", "SMILES", "Formula"])  # Write header
            for cid, mw, smiles, formula in candidates:
                writer.writerow([cid, mw, smiles, formula])  # Write data
        f.close()


def search_candidates_by_mw(input_data, pubchem_dict, candidates_folder, ppm_threshold):
    # 解析 CID 和 目标分子量
    target_mws = []
    for _, row in input_data.iterrows():
        cid = int(row[0])  # 第一列是物质 ID
        mw = float(row[1])  # 第二列是分子量
        target_mws.append((cid, mw))

    pubchem_mws = [
        (cid, entry["mw"], entry["smiles"], entry["formula"])
        for cid, entry in tqdm(pubchem_dict.items(), desc="Unpacking pubchem dict")
    ]
    print("len of pubchem_mws:", len(pubchem_mws))

    result = find_candidates_by_ppm_smiles(target_mws, pubchem_mws, ppm_threshold)

    save_results_to_csv(result, candidates_folder)

    print("All results saved!")


def filter_candidates_by_formula(
    input_folder, output_folder, is_filter_element=True, element_set=None
):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if is_filter_element:
        # 允许的元素集合（注意顺序处理避免误识别Cl为C和l）
        allowed_elements = set(element_set.split(","))
        # allowed_elements = {'C', 'H', 'O', 'N', 'P', 'S', 'F', 'Cl', 'Br', 'I'}

        # 用于正则提取分子式中的元素符号（如 C, Cl, Br）
        element_pattern = re.compile(r"Cl|Br|[A-Z][a-z]?")

        # 获取所有csv文件名
        csv_files = [f for f in os.listdir(input_folder) if f.endswith(".csv")]

        # 初始化进度条
        with tqdm(total=len(csv_files), desc="Filtering formulas") as pbar:
            for file in csv_files:
                # 读取文件
                file_path = os.path.join(input_folder, file)
                df = pd.read_csv(file_path)

                # 检查是否有Formula列
                if "Formula" not in df.columns:
                    tqdm.write(f"Skipped {file}: No 'Formula' column found.")
                    pbar.update(1)
                    continue

                # 过滤非法分子式
                def is_valid_formula(formula):
                    elements = element_pattern.findall(str(formula))
                    return all(e in allowed_elements for e in elements)

                df_filtered = df[df["Formula"].apply(is_valid_formula)]

                # 写入输出文件
                output_path = os.path.join(output_folder, file)
                df_filtered.to_csv(output_path, index=False)

                pbar.update(1)


def filter_candidates_by_smiles(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # 创建输出目录（如不存在）
    os.makedirs(output_folder, exist_ok=True)

    # 获取所有csv文件名
    csv_files = [f for f in os.listdir(input_folder) if f.endswith(".csv")]

    # 初始化进度条
    with tqdm(total=len(csv_files), desc="Filtering candidates smiles") as pbar:
        for file in csv_files:
            file_path = os.path.join(input_folder, file)
            try:
                df = pd.read_csv(file_path)
                if "SMILES" not in df.columns:
                    print(f"⚠️ Skipped {file}: No 'SMILES' column found.")
                    pbar.update(1)
                    continue

                # 保留不含 '.', '+', '-' 的行
                df_filtered = df[
                    ~df["SMILES"].fillna("").astype(str).str.contains(r"[.\+\-]")
                ]

                # 保存过滤后的数据
                output_path = os.path.join(output_folder, file)
                df_filtered.to_csv(output_path, index=False)

            except Exception as e:
                # print(f"❌ Error processing {file}: {e}")
                print(f"{file} is empty or has no valid data, skipped.")
                continue
            pbar.update(1)


def remove_empty_csv_files(folder_path):
    cnt = 0
    for filename in os.listdir(folder_path):
        if filename.endswith(".csv"):
            file_path = os.path.join(folder_path, filename)
            df = pd.read_csv(file_path)
            if df.empty:
                cnt += 1
                os.remove(file_path)
    print(f"Total empty files removed: {cnt}")


if __name__ == "__main__":
    start_time = time.time()
    args = parse_args()

    tmp_result_path = "tmp/"
    if not os.path.exists(tmp_result_path):
        os.makedirs(tmp_result_path)

    # 生成唯一 ID 和分子量的 CSV 文件
    molecular_network_df = pd.read_csv(args.edited_molecular_network_file)
    unique_ids_df = generate_ID_MW_file(
        molecular_network_df, os.path.join(tmp_result_path, "unique_ID_MW.csv")
    )

    # 加载 PubChem 数据库
    print("Loading pubchem database...")
    with open(args.pubchem_database_path, "rb") as f:
        pubchem_dict = pk.load(f)
        f.close()

    # 根据分子质量的ppm误差范围搜索 PubChem 数据库中符合条件的分子
    tmp_candidates_folder = tmp_result_path + "candidates"
    if not os.path.exists(tmp_candidates_folder):
        os.makedirs(tmp_candidates_folder)
    search_candidates_by_mw(
        unique_ids_df, pubchem_dict, tmp_candidates_folder, args.ppm_threshold
    )

    # 依据分子式进行元素过滤，得到candidates
    tmp_filtered_candidates_by_formula_folder = (
        tmp_result_path + "candidates_filtered_by_formula"
    )
    filter_candidates_by_formula(
        tmp_candidates_folder,
        tmp_filtered_candidates_by_formula_folder,
        args.is_filter_element,
        args.element_set,
    )

    # 依据SMILES进行过滤，得到candidates
    filter_candidates_by_smiles(
        tmp_filtered_candidates_by_formula_folder, args.candidates_folder
    )

    remove_empty_csv_files(args.candidates_folder)
    print("Spend time: ", time.time() - start_time)
