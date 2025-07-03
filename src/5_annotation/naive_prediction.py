import os
import time
import shutil
import argparse
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, DataStructs


def parse_args():
    parser = argparse.ArgumentParser(description="Navie prediction.")
    parser.add_argument(
        "--molecular_network_file",
        default="source_target_cor_edit.csv",
        type=str,
        required=True,
        help="input molecular network file containing source, target, corrrelation",
    )
    parser.add_argument(
        "--seednode_file",
        default="seednode.csv",
        type=str,
        required=True,
        help="the original seednode file containing ID and SMILES",
    )
    parser.add_argument(
        "--threshold_tanimoto_similarity",
        default=0.5,
        type=float,
        required=True,
        help="the threshold for tanimoto similarity to filter the results",
    )
    return parser.parse_args()


def pick_cycle_annotation(
    unique_folder,
    not_unique_folder,
    naive_prediction_folder,
    threshold_tanimoto_similarity,
):
    # 初始化结果列表
    result_rows = []

    # 获取所有CSV文件
    unique_files = [f for f in os.listdir(unique_folder) if f.endswith(".csv")]
    not_unique_files = [f for f in os.listdir(not_unique_folder) if f.endswith(".csv")]

    # 将unique_files和not_unique_files copy到naive_prediction_folder
    for file in unique_files + not_unique_files:
        src_path = os.path.join(
            unique_folder if file in unique_files else not_unique_folder, file
        )
        dest_path = os.path.join(naive_prediction_folder, file)
        shutil.copy(src_path, dest_path)

    # 遍历unique_files
    for file in unique_files:
        file_path = os.path.join(naive_prediction_folder, file)
        df = pd.read_csv(file_path)
        # 从file中提取ID
        file_id = file.split(".csv")[0]
        # 从df中找出score值最大的对应的ID
        if "score" in df.columns:
            max_score = df["score"].max()
            smiles = df.loc[df["score"] == max_score, "SMILES"].values[0]
            if max_score < threshold_tanimoto_similarity:
                continue
        else:
            continue
        result_rows.append({"ID": file_id, "SMILES": smiles})

    # 遍历not_unique_files
    for file in not_unique_files:
        file_path = os.path.join(naive_prediction_folder, file)
        df = pd.read_csv(file_path)
        # 从file中提取ID
        file_id = file.split(".csv")[0]
        # 从df中找出score值最大的对应的ID
        if "weighted_score" in df.columns:
            max_score = df["weighted_score"].max()
            smiles = df.loc[df["weighted_score"] == max_score, "SMILES"].values[0]
            if max_score < threshold_tanimoto_similarity:
                continue
        else:
            continue
        result_rows.append({"ID": file_id, "SMILES": smiles})

    # 将结果转换为DataFrame
    final_df = pd.DataFrame(result_rows)
    return final_df


if __name__ == "__main__":
    args = parse_args()
    seednode_file = args.seednode_file
    threshold_tanimoto_similarity = args.threshold_tanimoto_similarity

    # 读取 all_nodes 集合
    molecular_network_df = pd.read_csv(args.molecular_network_file)
    source_nodes = molecular_network_df["source"].tolist()
    target_nodes = molecular_network_df["target"].tolist()
    all_nodes = set(source_nodes + target_nodes)

    unique_folder = (
        "tmp/Seednode_and_Targetnode_Morgan_Similarity_score_split_unique_Top10"
    )
    not_unique_folder = (
        "tmp/Seednode_and_Targetnode_Morgan_Similarity_score_split_not_unique_Top10"
    )
    naive_prediction_folder = "tmp/naive_prediction_results"
    os.makedirs(naive_prediction_folder, exist_ok=True)

    # 运行 naive prediction
    final_df = pick_cycle_annotation(
        unique_folder,
        not_unique_folder,
        naive_prediction_folder,
        threshold_tanimoto_similarity,
    )

    # 生成下一轮的seednode file
    seednode_df = pd.read_csv(seednode_file)
    seednode_df = seednode_df[["ID", "SMILES"]]
    merged_df = pd.concat([seednode_df, final_df], ignore_index=True)

    merged_df.to_csv("tmp/naive_cycle_seednode.csv", index=False)
