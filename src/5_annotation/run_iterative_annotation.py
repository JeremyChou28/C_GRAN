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
        "--num_containers",
        default=10,
        type=int,
        required=True,
        help="the number of Docker containers to run in parallel",
    )
    parser.add_argument(
        "--tolerance",
        default=0.1,
        type=float,
        required=True,
        help="the tolerance for the modified cosine similarity calculation",
    )
    parser.add_argument(
        "--energy_level",
        default=0,
        type=int,
        required=True,
        choices=[0, 1, 2],
        help="the energy level of the spectrum data predicted by CFM-ID, 0 represents low, 1 represents medium, and 2 represents high energy level",
    )
    parser.add_argument(
        "--ion_mode",
        default="positive",
        type=str,
        required=True,
        choices=["positive", "negative"],
        help="the ion mode of the spectrum data predicted by CFM-ID, positive or negative",
    )
    parser.add_argument(
        "--modified_cosine_similarity_threshold",
        default=0.5,
        type=float,
        required=True,
        help="the threshold for modified cosine similarity to filter the results",
    )
    parser.add_argument(
        "--spectrum_file",
        default="./test_files/compounds_spectrum.mgf",
        type=str,
        required=True,
        help="the mgf file containing target node spectra",
    )
    parser.add_argument(
        "--top_k",
        type=int,
        required=True,
        default=10,
        help="the top k candidates for annotation",
    )
    parser.add_argument(
        "--max_iterations",
        default=100,
        type=int,
        required=True,
        help="the maximum number of iterations to run the annotation process",
    )
    return parser.parse_args()


def merge_cfmid_results(cfmid_score_results_path, output_file, cfmid_seednodes_path):
    id_to_round = {}
    round_num = int(
        file.replace("annotation_with_cfmid_seednode_round", "").replace(".csv", "")
    )
    round_df = pd.read_csv(os.path.join(cfmid_seednodes_path, file))
    for id_val in round_df["ID"].dropna().unique():
        id_to_round[str(id_val)] = round_num  # 保证 key 是字符串格式

    all_files = [f for f in os.listdir(cfmid_score_results_path) if f.endswith(".csv")]
    df_list = []

    for file in all_files:
        file_path = os.path.join(cfmid_score_results_path, file)
        df = pd.read_csv(file_path)
        # 提取ID并添加为第一列
        id_value = file.split(".")[0]  # 假设文件名格式为 ID.csv
        # 判断id_value是否存在于 id_to_round 中
        if id_value not in id_to_round:
            continue

        df.insert(0, "ID", id_value)  # 插入到第一列

        df.rename(columns={"Seednode": "Seed Node"}, inplace=True)
        df.rename(columns={"CFM-ID_score": "Score"}, inplace=True)
        df.rename(columns={"MW": "MonoIsotopic Weight"}, inplace=True)

        columns_to_keep = [
            "ID",
            "Seed Node",
            "CID",
            "MonoIsotopic Weight",
            "SMILES",
            "Formula",
            "Score",
        ]
        df_final = df[columns_to_keep].copy()

        df_final.sort_values(by="Score", ascending=False, inplace=True)

        # 添加 Round 列作为第一列
        round_value = id_to_round.get(str(id_value), "")
        df_final.insert(0, "Round", round_value)

        df_list.append(df_final)

    # 合并所有DataFrame
    final_df = pd.concat(df_list, ignore_index=True)
    final_df.to_csv(output_file, index=False)
    print(f"All results merged into {output_file}")


if __name__ == "__main__":
    start_time = time.time()
    args = parse_args()
    # 初始的seednode file
    seednode_file = args.seednode_file
    if not os.path.exists("tmp"):
        os.makedirs("tmp")
    shutil.copy(seednode_file, "tmp/annotation_with_cfmid_seednode_round0.csv")

    last_cycle_seednode_df = pd.read_csv(seednode_file)
    last_cycle_ids = set(last_cycle_seednode_df["ID"].astype(int).tolist())

    # 读取 all_nodes 集合
    molecular_network_df = pd.read_csv(args.edited_molecular_network_file)
    source_nodes = molecular_network_df["source"].tolist()
    target_nodes = molecular_network_df["target"].tolist()

    # 循环直到 seednode 中包含所有节点
    max_rounds = args.max_iterations  # 防止死循环，你也可以去掉
    round_num = 0

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
                "cfmid_prediction.py",
                "--num_containers",
                str(args.num_containers),
                "--tolerance",
                str(args.tolerance),
                "--energy_level",
                str(args.energy_level),
                "--ion_mode",
                args.ion_mode,
                "--spectrum_file",
                args.spectrum_file,
                "--top_k",
                str(args.top_k),
                "--modified_cosine_similarity_threshold",
                str(args.modified_cosine_similarity_threshold),
                "--round_num",
                str(round_num),
            ],
            check=True,
        )
        subprocess.run(
            [
                "python",
                "postprocess.py",
                "--seednode_file",
                seednode_file,
                "--round_num",
                str(round_num),
            ],
            check=True,
        )

        seednode_file = f"tmp/annotation_with_cfmid_seednode_round{round_num}.csv"
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
    merge_cfmid_results(
        cfmid_score_results_path="tmp/cfmid_score_results/",
        output_file="final_cfmid_annotation_results.csv",
        cfmid_seednodes_path=f"tmp/annotation_with_cfmid_seednode_round{round_num}.csv",
    )
    print("Spend time: ", time.time() - start_time)
