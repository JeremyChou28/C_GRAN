import os
import argparse
import subprocess
import pandas as pd

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
        "--threshold_tanimoto_similarity",
        default=0.5,
        type=float,
        required=True,
        help="the threshold for tanimoto similarity to filter the results",
    )
    return parser.parse_args()


def merge_naive_results(tmp_folder, output_file):
    all_files = [f for f in os.listdir(tmp_folder) if f.endswith('.csv')]
    df_list = []

    for file in all_files:
        file_path = os.path.join(tmp_folder, file)
        df = pd.read_csv(file_path)

        # 提取 ID（文件名去掉 .csv）
        id_value = os.path.splitext(file)[0]

        # 统一字段：确保包含 weighted_score 列（若无则补 NaN）
        if 'weighted_score' not in df.columns:
            df['weighted_score'] = pd.NA  # 或 None

        # 保留并重排你需要的列
        df = df[['CID', 'MW', 'SMILES', 'Formula', 'score', 'weighted_score']]
        df.insert(0, 'ID', id_value)  # 将 ID 插入第一列

        df_list.append(df)

    # 合并所有DataFrame
    final_df = pd.concat(df_list, ignore_index=True)
    final_df.to_csv(output_file, index=False)
    print(f"All results merged into {output_file}")

if __name__ == "__main__":
    args=parse_args()
    # 初始的seednode file
    seednode_file = args.seednode_file
    
    # 读取 all_nodes 集合
    molecular_network_df = pd.read_csv(args.molecular_network_file)
    source_nodes = molecular_network_df['source'].tolist()
    target_nodes = molecular_network_df['target'].tolist()
    all_nodes = set(source_nodes + target_nodes)

    # 循环直到 seednode 中包含所有节点
    max_rounds = 100  # 防止死循环，你也可以去掉
    round_num = 0

    while round_num < max_rounds:
        round_num += 1
        print(f"\n🔁 Round {round_num} running...")

        # 依次运行三个脚本
        subprocess.run(["python", "preprocess.py", "--molecular_network_file", args.molecular_network_file, "--seednode_file", seednode_file, "--candidates_folder", args.candidates_folder], check=True)
        subprocess.run(["python", "naive_prediction.py", "--molecular_network_file", args.molecular_network_file, "--seednode_file", seednode_file,"--threshold_tanimoto_similarity", str(args.threshold_tanimoto_similarity)], check=True)

        seednode_file='tmp/naive_cycle_seednode.csv'
        # 读取 seednode 文件，检查 ID 集合
        try:
            seednode_df = pd.read_csv(seednode_file)
            current_ids = set(seednode_df['ID'].astype(int).tolist())
        except Exception as e:
            print(f"❌ Failed to read seednode.csv: {e}")
            break

        # 判断是否完成
        if current_ids == all_nodes:
            print("🎉 注释完成！")
            break
    else:
        print("⚠️ 达到最大循环次数，仍未完成。")

    # 输出最终的注释结果
    merge_naive_results(tmp_folder='tmp/naive_prediction_results', output_file='final_naive_annotation_results.csv')