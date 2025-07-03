import os
import time
import shutil
import argparse
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, DataStructs


def parse_args():
    parser = argparse.ArgumentParser(description="Molecular Networking Annotation.")
    parser.add_argument(
        "--molecular_network_file",
        default="source_target_cor_edit.csv",
        type=str,
        required=True,
        help="input molecular network file containing source, target, corrrelation",
    )
    parser.add_argument(
        "--seednode_file",
        type=str,
        required=True,
        help="seed node file containing ID and SMILES",
    )
    parser.add_argument(
        "--candidates_folder",
        type=str,
        help="the candidates folder to save the search results",
    )
    parser.add_argument(
        "--top_k",
        type=int,
        required=True,
        default=10,
        help="the top k candidates for annotation",
    )
    return parser.parse_args()


def generate_seednode_and_targetnode(
    molecular_network_df, seednode_df, output_filename
):
    # 提取种子物质 ID 集合
    seed_ids = set(seednode_df["ID"].astype(int))  # 确保为整数类型以匹配

    # 初始化列表用于存储结果
    results = []

    # 使用 tqdm 包装 DataFrame 行迭代器
    for _, row in tqdm(
        molecular_network_df.iterrows(),
        total=len(molecular_network_df),
        desc="Filtering connected nodes",
    ):
        source = int(row["source"])
        target = int(row["target"])
        corr = row["Correlation"]  # 提取 Correlation 值

        source_in_seed = source in seed_ids
        target_in_seed = target in seed_ids

        # 仅保留一端在种子集合中，另一端不在
        if source_in_seed and not target_in_seed:
            results.append(
                {"Seednode": source, "Targetnode": target, "Correlation": corr}
            )
        elif target_in_seed and not source_in_seed:
            results.append(
                {"Seednode": target, "Targetnode": source, "Correlation": corr}
            )

    # 转换为 DataFrame 并按 Seednode 升序排序
    result_df = pd.DataFrame(results)
    result_df.sort_values(by="Seednode", inplace=True)

    # 保存输出文件
    result_df.to_csv(output_filename, index=False)
    return result_df


def calculate_tanimoto_similarity(
    seednode_df, edges_df, candidates_folder, output_folder
):
    # 读取 Seednode 的 SMILES 信息
    seednode_smiles_dict = dict(zip(seednode_df["ID"], seednode_df["SMILES"]))

    def tanimoto_similarity(smiles1, smiles2):
        """计算两个 SMILES 之间的 Tanimoto 相似性"""
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        if mol1 is None or mol2 is None:
            return None  # 处理 SMILES 解析失败的情况
        fp1 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=2048)
        fp2 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=2048)
        return DataStructs.TanimotoSimilarity(fp1, fp2)

    # 计算并保存相似性分数
    for _, row in tqdm(
        edges_df.iterrows(), total=len(edges_df), desc="Processing Edges"
    ):
        target_id, seed_id = int(row["Targetnode"]), int(row["Seednode"])

        # 获取 Seednode 的 SMILES
        seed_smiles = seednode_smiles_dict.get(seed_id)
        if not seed_smiles:
            continue  # 如果 Seednode 没有 SMILES，跳过

        # 读取 Targetnode 文件
        target_file = os.path.join(candidates_folder, f"{target_id}.csv")
        if not os.path.exists(target_file):
            continue  # 如果没有该 Targetnode 文件，跳过

        df_target = pd.read_csv(target_file)
        if "SMILES" not in df_target.columns:
            continue  # 确保文件中包含 SMILES 列

        # 计算相似性并添加到 DataFrame
        df_target["score"] = df_target["SMILES"].apply(
            lambda x: tanimoto_similarity(x, seed_smiles)
        )

        # 按 score 降序排序
        df_target = df_target.sort_values(by="score", ascending=False)

        # 保存到新文件
        output_filename = os.path.join(output_folder, f"{target_id}_{seed_id}.csv")
        df_target.to_csv(output_filename, index=False)

    print("相似性计算完成，结果已保存到 output 目录。")


def group_target_nodes(seed_target_df, output_unique, output_repeated):
    # 使用 tqdm 展示拆分进度
    with tqdm(total=1, desc="Splitting based on Targetnode count") as pbar:
        # 计算每个 Targetnode 出现次数
        target_counts = seed_target_df["Targetnode"].value_counts()

        # 获取只出现一次和多次的 Targetnode ID 集合
        unique_targets = set(target_counts[target_counts == 1].index)
        repeated_targets = set(target_counts[target_counts > 1].index)

        # 拆分为两个 DataFrame
        df_unique = seed_target_df[
            seed_target_df["Targetnode"].isin(unique_targets)
        ].copy()
        df_repeated = seed_target_df[
            seed_target_df["Targetnode"].isin(repeated_targets)
        ].copy()

        pbar.update(1)

    # 保存两个文件
    df_unique.to_csv(output_unique, index=False)
    df_repeated.to_csv(output_repeated, index=False)


def split_target_node(
    unique_csv,
    not_unique_csv,
    seed_target_tanimoto_folder,
    unique_folder,
    not_unique_folder,
):

    # 读取两个csv文件
    df_unique = pd.read_csv(unique_csv)
    df_not_unique = pd.read_csv(not_unique_csv)

    # 构建命名集合（文件名：'Targetnode_Seednode.csv'）
    def build_filename_set(df):
        return set(
            f"{int(row['Targetnode'])}_{int(row['Seednode'])}.csv"
            for _, row in df.iterrows()
        )

    unique_files = build_filename_set(df_unique)
    not_unique_files = build_filename_set(df_not_unique)

    # 遍历原始文件夹中的所有csv文件，按匹配规则分类拷贝
    all_files = [
        f for f in os.listdir(seed_target_tanimoto_folder) if f.endswith(".csv")
    ]

    with tqdm(total=len(all_files), desc="Classifying CSV files") as pbar:
        for file in all_files:
            source_path = os.path.join(seed_target_tanimoto_folder, file)

            if file in unique_files:
                shutil.copy(source_path, os.path.join(unique_folder, file))
            elif file in not_unique_files:
                shutil.copy(source_path, os.path.join(not_unique_folder, file))
            # 忽略未匹配文件
            pbar.update(1)


def pickle_topk_unique_seednode(unique_folder, topk):
    output_folder = f"tmp/Seednode_and_Targetnode_Morgan_Similarity_score_split_unique_Top{topk}"  # 输出文件夹

    # 创建输出文件夹
    os.makedirs(output_folder, exist_ok=True)

    # 获取所有以 "_xxx.csv" 命名的csv文件
    csv_files = [
        f for f in os.listdir(unique_folder) if f.endswith(".csv") and "_" in f
    ]

    # 初始化进度条
    with tqdm(total=len(csv_files), desc="Filtering top score rows") as pbar:
        for filename in csv_files:
            filepath = os.path.join(unique_folder, filename)
            try:
                df = pd.read_csv(filepath)
                if "score" not in df.columns:
                    tqdm.write(f"Skipped {filename}: no 'score' column")
                    pbar.update(1)
                    continue

                # 获取前k名及并列的 score
                scores = df["score"].tolist()
                if len(scores) <= topk:
                    df_filtered = df
                else:
                    threshold_score = sorted(scores, reverse=True)[topk - 1]
                    df_filtered = df[df["score"] >= threshold_score]

                # 输出文件名：取第一个数值部分整数形式
                prefix = filename.split("_")[0]  # 例如 '262.0'
                output_filename = f"{int(float(prefix))}.csv"
                output_path = os.path.join(output_folder, output_filename)

                # 保存（覆盖）
                df_filtered.to_csv(output_path, index=False)

            except Exception as e:
                tqdm.write(f"Error processing {filename}: {str(e)}")
            pbar.update(1)


def pickle_topk_not_unique(not_unique_csv, not_unique_folder, topk):
    # 路径设置（请修改为你实际的路径）
    output_folder = f"tmp/Seednode_and_Targetnode_Morgan_Similarity_score_split_not_unique_Top{topk}"  # 输出文件夹

    os.makedirs(output_folder, exist_ok=True)

    # 读取文件A，建立 (Targetnode, Seednode) → Correlation 映射
    df_corr = pd.read_csv(not_unique_csv)

    # 构建字典索引
    correlation_dict = {
        (int(row["Targetnode"]), int(row["Seednode"])): row["Correlation"]
        for _, row in df_corr.iterrows()
    }

    # 构建：Targetnode → 多文件 DataFrame 聚合表
    targetnode_groups = {}

    # 遍历文件夹B
    csv_files = [
        f for f in os.listdir(not_unique_folder) if f.endswith(".csv") and "_" in f
    ]

    with tqdm(total=len(csv_files), desc="Processing files") as pbar:
        for filename in csv_files:
            try:
                # 解析文件名：Targetnode_Seednode.csv
                base_name = filename[:-4]  # 去掉 .csv
                target_str, seed_str = base_name.split("_")
                targetnode = eval(target_str)
                seednode = eval(seed_str)

                # 读取文件并添加 Correlation 加权列
                file_path = os.path.join(not_unique_folder, filename)
                df = pd.read_csv(file_path)
                if "score" not in df.columns:
                    tqdm.write(f"Skipped {filename}: no 'score' column.")
                    pbar.update(1)
                    continue

                # 获取对应 Correlation 值
                correlation = correlation_dict.get((targetnode, seednode))
                if correlation is None:
                    tqdm.write(f"Skipped {filename}: no Correlation match.")
                    pbar.update(1)
                    continue

                # 加权打分
                df["weighted_score"] = df["score"] * correlation
                df["Targetnode"] = targetnode
                df["Seednode"] = seednode

                # 收集该 Targetnode 的所有行
                if targetnode not in targetnode_groups:
                    targetnode_groups[targetnode] = []
                targetnode_groups[targetnode].append(df)

            except Exception as e:
                tqdm.write(f"Error processing {filename}: {str(e)}")
            pbar.update(1)

    # 对每个 Targetnode 聚合并筛选 top-N weighted_score
    for targetnode, df_list in tqdm(
        targetnode_groups.items(), desc="Filtering and saving"
    ):
        combined_df = pd.concat(df_list, ignore_index=True)

        # 按 weighted_score 降序排列
        combined_df.sort_values(by="weighted_score", ascending=False, inplace=True)
        top_scores = combined_df["weighted_score"].tolist()

        # 处理并列：找前k名及相同分数
        if len(top_scores) <= topk:
            df_top = combined_df
        else:
            threshold = top_scores[topk - 1]
            df_top = combined_df[combined_df["weighted_score"] >= threshold]

        # 输出文件名：Targetnode 整数位.csv
        output_name = f"{int(float(targetnode))}.csv"
        output_path = os.path.join(output_folder, output_name)
        df_top.to_csv(output_path, index=False)


if __name__ == "__main__":
    start_time = time.time()
    args = parse_args()

    tmp_result_path = "tmp/"
    if not os.path.exists(tmp_result_path):
        os.makedirs(tmp_result_path)

    # 获取seed node和target node
    molecular_network_df = pd.read_csv(args.molecular_network_file)
    seednode_df = pd.read_csv(args.seednode_file)
    edges_df = generate_seednode_and_targetnode(
        molecular_network_df,
        seednode_df,
        os.path.join(tmp_result_path, "Seednode_and_Targetnode.csv"),
    )

    # 计算tanimoto similarity
    calculate_tanimoto_similarity(
        seednode_df,
        edges_df,
        args.candidates_folder,
        os.path.join(
            tmp_result_path, "Seednode_and_Targetnode_Morgan_Similarity_score"
        ),
    )

    # 检查种子节点
    seed_target_df = pd.read_csv(
        os.path.join(tmp_result_path, "Seednode_and_Targetnode.csv")
    )
    output_unique = tmp_result_path + "Seednode_and_Targetnode_unique_seednode.csv"
    output_repeated = (
        tmp_result_path + "Seednode_and_Targetnode_not_unique_seednodes.csv"
    )
    group_target_nodes(seed_target_df, output_unique, output_repeated)

    unique_csv = tmp_result_path + "Seednode_and_Targetnode_unique_seednode.csv"
    not_unique_csv = (
        tmp_result_path + "Seednode_and_Targetnode_not_unique_seednodes.csv"
    )
    seed_target_tanimoto_folder = (
        tmp_result_path + "Seednode_and_Targetnode_Morgan_Similarity_score"
    )
    unique_folder = "tmp/Seednode_and_Targetnode_Morgan_Similarity_score_split_unique"
    not_unique_folder = (
        "tmp/Seednode_and_Targetnode_Morgan_Similarity_score_split_not_unique"
    )
    os.makedirs(seed_target_tanimoto_folder, exist_ok=True)
    os.makedirs(unique_folder, exist_ok=True)
    os.makedirs(not_unique_folder, exist_ok=True)

    split_target_node(
        unique_csv,
        not_unique_csv,
        seed_target_tanimoto_folder,
        unique_folder,
        not_unique_folder,
    )

    pickle_topk_unique_seednode(unique_folder, args.top_k)
    pickle_topk_not_unique(not_unique_csv, not_unique_folder, args.top_k)

    print("Spend time: ", time.time() - start_time)
