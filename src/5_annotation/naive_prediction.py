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
        "--edited_molecular_network_file",
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
        "--tanimoto_similarity_threshold",
        default=0.5,
        type=float,
        required=True,
        help="the threshold for tanimoto similarity to filter the results",
    )
    parser.add_argument(
        "--top_k",
        type=int,
        required=True,
        default=10,
        help="the top k candidates for annotation",
    )
    parser.add_argument(
        "--round",
        type=int,
        default=0,
        help="the round of annotation, used to name the output file",
    )
    return parser.parse_args()


def pick_cycle_annotation(
    unique_folder,
    not_unique_folder,
    naive_prediction_folder,
    tanimoto_similarity_threshold,
):
    # 初始化结果列表
    result_rows = []

    # 获取所有CSV文件
    unique_files = [f for f in os.listdir(unique_folder) if f.endswith(".csv")]
    not_unique_files = [f for f in os.listdir(not_unique_folder) if f.endswith(".csv")]

    # 遍历unique_files
    for file in unique_files:
        file_path = os.path.join(unique_folder, file)
        df = pd.read_csv(file_path)
        # 从file中提取ID
        file_id = file.split(".csv")[0]
        # 从df中找出score值最大的对应的ID
        if "score" in df.columns:
            max_score = df["score"].max()
            smiles = df.loc[df["score"] == max_score, "SMILES"].values[0]
            if max_score < tanimoto_similarity_threshold:
                continue
        else:
            continue
        result_rows.append({"ID": file_id, "SMILES": smiles})

    # 遍历not_unique_files
    for file in not_unique_files:
        file_path = os.path.join(not_unique_folder, file)
        df = pd.read_csv(file_path)
        # 从file中提取ID
        file_id = file.split(".csv")[0]
        # 从df中找出score值最大的对应的ID
        if "weighted_score" in df.columns:
            max_score = df["weighted_score"].max()
            smiles = df.loc[df["weighted_score"] == max_score, "SMILES"].values[0]
            if max_score < tanimoto_similarity_threshold:
                continue
        else:
            continue
        result_rows.append({"ID": file_id, "SMILES": smiles})


    # 将结果转换为DataFrame
    final_df = pd.DataFrame(result_rows)

    # 将df中的ID列的文件copy到naive_prediction_folder
    for index in final_df['ID'].tolist():
        file_name = f"{index}.csv"
        src_path = os.path.join(unique_folder, file_name)
        if not os.path.exists(src_path):
            src_path = os.path.join(not_unique_folder, file_name)
        dest_path = os.path.join(naive_prediction_folder, file_name)
        shutil.copy(src_path, dest_path)

    return final_df


if __name__ == "__main__":
    args = parse_args()
    seednode_file = args.seednode_file
    tanimoto_similarity_threshold = args.tanimoto_similarity_threshold

    # 读取 all_nodes 集合
    molecular_network_df = pd.read_csv(args.edited_molecular_network_file)
    source_nodes = molecular_network_df["source"].tolist()
    target_nodes = molecular_network_df["target"].tolist()
    all_nodes = set(source_nodes + target_nodes)

    unique_folder = (
        f"tmp/Seednode_and_Targetnode_Morgan_Similarity_score_split_unique_Top{args.top_k}"
    )
    not_unique_folder = (
        f"tmp/Seednode_and_Targetnode_Morgan_Similarity_score_split_not_unique_Top{args.top_k}"
    )
    naive_prediction_folder = "tmp/naive_prediction_results"
    os.makedirs(naive_prediction_folder, exist_ok=True)

    # 运行 naive prediction
    final_df = pick_cycle_annotation(
        unique_folder,
        not_unique_folder,
        naive_prediction_folder,
        tanimoto_similarity_threshold,
    )

    # 生成下一轮的seednode file
    seednode_df = pd.read_csv(seednode_file)
    seednode_df = seednode_df[["ID", "SMILES"]]
    merged_df = pd.concat([seednode_df, final_df], ignore_index=True)

    merged_df.to_csv(f"tmp/naive_annotation_seednode_round{args.round}.csv", index=False)
