import os
import time
import argparse
import subprocess
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Molecular Networking Annotation.")
    parser.add_argument(
        "--edited_molecular_network_file",
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
        "--tanimoto_similarity_threshold",
        default=0.5,
        type=float,
        required=True,
        help="the threshold for tanimoto similarity to filter the results",
    )
    parser.add_argument(
        "--max_iterations",
        default=100,
        type=int,
        required=True,
        help="the maximum number of iterations to run the annotation process",
    )
    parser.add_argument(
        "--top_k",
        type=int,
        required=True,
        default=10,
        help="the top k candidates for annotation",
    )
    return parser.parse_args()


def merge_naive_results(
    naive_prediction_results_path, output_file, naive_seednodes_path, iterations
):
    all_files = [
        f for f in os.listdir(naive_prediction_results_path) if f.endswith(".csv")
    ]
    df_list = []

    # Step 1: 构建 ID -> Round 映射
    id_to_round = {}
    for round_num in range(1, iterations + 1):
        file = f"naive_annotation_seednode_round{round_num}.csv"
        file_path = os.path.join(naive_seednodes_path, file)
        df = pd.read_csv(file_path)

        for node_id in df["ID"]:
            if str(node_id) not in id_to_round:
                id_to_round[str(node_id)] = round_num

    # Step 2: 合并每个结果文件
    for file in all_files:
        file_path = os.path.join(naive_prediction_results_path, file)
        df = pd.read_csv(file_path)

        # 当前文件的 ID（由文件名去除 .csv 得到）
        id_value = os.path.splitext(file)[0]

        # 标准化列名
        df.rename(columns={"Seednode": "Seed Node"}, inplace=True)
        df.rename(columns={"MW": "MonoIsotopic Weight"}, inplace=True)

        # 标准化 Score 列
        if "weighted_score" in df.columns:
            df["score"] = df["weighted_score"]
        df.rename(columns={"score": "Score"}, inplace=True)

        # 保留并排列需要的列
        df = df[
            ["Seed Node", "CID", "MonoIsotopic Weight", "SMILES", "Formula", "Score"]
        ]
        df.insert(0, "ID", id_value)  # 插入 ID 列
        df["CID"] = "https://pubchem.ncbi.nlm.nih.gov/compound/" + df["CID"].astype(str)
        df["Seed Node"] = df["Seed Node"].apply(
            lambda x: "" if pd.isna(x) else str(int(x))
        )

        # 添加 Round 列作为第一列
        round_value = id_to_round.get(id_value, "")
        df.insert(0, "Round", [round_value] * len(df))

        df_list.append(df)

    # Step 3: 合并所有 DataFrame
    final_df = pd.concat(df_list, ignore_index=True)
    final_df.to_csv(output_file, index=False)
    print(f"✅ All results merged into {output_file}")


if __name__ == "__main__":
    start_time = time.time()
    args = parse_args()
    # 初始的seednode file
    seednode_file = args.seednode_file

    last_cycle_seednode_df = pd.read_csv(seednode_file)
    last_cycle_ids = set(last_cycle_seednode_df["ID"].astype(int).tolist())

    # 读取 all_nodes 集合
    molecular_network_df = pd.read_csv(args.edited_molecular_network_file)
    source_nodes = molecular_network_df["source"].tolist()
    target_nodes = molecular_network_df["target"].tolist()

    # 循环直到 seednode 中包含所有节点
    max_rounds = args.max_iterations  # 防止死循环，你也可以去掉

    while round_num < max_rounds:
        round_time = time.time()
        round_num += 1
        print(f"\n🔁 Round {round_num} running...")

        # 依次运行三个脚本
        subprocess.run(
            [
                "python",
                "preprocess.py",
                "--edited_molecular_network_file",
                args.edited_molecular_network_file,
                "--seednode_file",
                seednode_file,
                "--candidates_folder",
                args.candidates_folder,
                "--top_k",
                str(args.top_k),
                "--round_num",
                str(round_num),
            ],
            check=True,
        )
        subprocess.run(
            [
                "python",
                "naive_prediction.py",
                "--edited_molecular_network_file",
                args.edited_molecular_network_file,
                "--seednode_file",
                seednode_file,
                "--tanimoto_similarity_threshold",
                str(args.tanimoto_similarity_threshold),
                "--top_k",
                str(args.top_k),
                "--round_num",
                str(round_num),
            ],
            check=True,
        )

        seednode_file = f"tmp/naive_annotation_seednode_round{round_num}.csv"
        # 读取 seednode 文件，检查 ID 集合
        try:
            seednode_df = pd.read_csv(seednode_file)
            current_ids = set(seednode_df["ID"].astype(int).tolist())
        except Exception as e:
            print(f"❌ Failed to read seednode.csv: {e}")
            break

        # 判断是否完成
        if last_cycle_ids < current_ids:
            last_cycle_ids = current_ids
        else:
            round_num -= 1  # 如果没有新增ID，则注释完成，则round_num减1，确保最终输出的迭代次数是实际的
            print("🎉 注释完成！")
            break
        print("Round time: ", time.time() - round_time)
    else:
        print("⚠️ 达到最大循环次数，仍未完成。")

    # 输出最终的注释结果
    merge_naive_results(
        naive_prediction_results_path="tmp/naive_prediction_results/",
        output_file="final_naive_annotation_results.csv",
        naive_seednodes_path=f"tmp/",
        iterations=round_num,
    )
    print("Spend time: ", time.time() - start_time)
