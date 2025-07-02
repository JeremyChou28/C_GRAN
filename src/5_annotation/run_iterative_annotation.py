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
        "--threshold_modified_cosine_similarity",
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
        "--max_iterations",
        default=100,
        type=int,
        required=True,
        help="the maximum number of iterations to run the annotation process",
    )
    return parser.parse_args()


def merge_cfmid_results(tmp_folder, output_file):
    all_files = [f for f in os.listdir(tmp_folder) if f.endswith('.csv')]
    df_list = []

    for file in all_files:
        file_path = os.path.join(tmp_folder, file)
        df = pd.read_csv(file_path)
        # 提取ID并添加为第一列
        id_value = file.split('.')[0]  # 假设文件名格式为 ID.csv
        df.insert(0, 'ID', id_value)  # 插入到第一列
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
    max_rounds = args.max_iterations  # 防止死循环，你也可以去掉
    round_num = 0

    while round_num < max_rounds:
        round_num += 1
        print(f"\n🔁 Round {round_num} running...")

        # 依次运行三个脚本
        subprocess.run(["python", "preprocess.py", "--molecular_network_file", args.molecular_network_file, "--seednode_file", seednode_file, "--candidates_folder", args.candidates_folder], check=True)
        subprocess.run(["python", "cfmid_prediction.py", "--num_containers", str(args.num_containers), "--tolerance", str(args.tolerance),"--energy_level", str(args.energy_level), "--spectrum_file", args.spectrum_file], check=True)
        subprocess.run(["python", "postprocess.py","--seednode_file",seednode_file,"--threshold_modified_cosine_similarity",str(args.threshold_modified_cosine_similarity)], check=True)

        seednode_file='tmp/seednode.csv'
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
    merge_cfmid_results(tmp_folder='tmp/cfmid_score_results', output_file='final_annotation_results.csv')