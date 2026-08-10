import os
import time
import pandas as pd
from tqdm import tqdm
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter high correlation compounds and generate high correlated nodes."
    )
    # online parameters
    parser.add_argument(
        "--seednode_file",
        default="seednode.csv",
        type=str,
        required=True,
        help="input seednode file",
    )
    parser.add_argument(
        "--correlation_threshold",
        type=float,
        default=None,
        help="Correlation threshold for filtering compounds. If not set, no correlation filtering is applied.",
    )
    # offline parameters
    parser.add_argument(
        "--correlation_file",
        type=str,
        required=True,
        help="input correlation file",
    )
    parser.add_argument(
        "--tmp_name",
        default="test",
        type=str,
        required=True,
        help="directory to save temporary files",
    )
    return parser.parse_args()


def filter_according_to_substance(correlation_df, seednode_df, output_filename):
    # 提取第二个CSV文件中的ID列
    substances = seednode_df["ID"].tolist()

    # 筛选Substance 1或者Substance 2在ID中的行
    filtered_df = correlation_df[
        (correlation_df["Substance 1"].isin(substances))
        | (correlation_df["Substance 2"].isin(substances))
    ].copy()

    # 
    filtered_df["Substance 1"] = (
        filtered_df["Substance 1"].astype("string").str.strip()
    )
    filtered_df["Substance 2"] = (
        filtered_df["Substance 2"].astype("string").str.strip()
    )

    # 保存筛选后的结果到新的CSV文件，添加进度条
    with tqdm(total=filtered_df.shape[0]) as pbar:
        filtered_df.to_csv(output_filename, index=False)
        pbar.update(filtered_df.shape[0])
    return filtered_df


def generate_high_correlated_nodes(df, output_filename, threshold):
    # 确保 Correlation 是数值类型（遇到脏数据转成 NaN）
    df = df.copy()
    df["Correlation"] = pd.to_numeric(df["Correlation"], errors="coerce")

    # 如果用户没有设置阈值：不做相关性过滤，直接用全部 df
    if threshold is None:
        filtered_df = df
    else:
        # threshold 是数值时才过滤；同时把 NaN 的相关性自动排除
        filtered_df = df[df["Correlation"].notna() & (df["Correlation"] > float(threshold))]

    # 提取 Substance 1/2 的所有 ID
    with tqdm(total=2, desc="Processing IDs") as pbar:
        substances = pd.concat([filtered_df["Substance 1"], filtered_df["Substance 2"]], ignore_index=True)
        pbar.update(1)

        # 去重、排序
        unique_ids = pd.DataFrame(sorted(pd.unique(substances.dropna())), columns=["ID"])
        pbar.update(1)

    unique_ids.to_csv(output_filename, index=False)
    print(f"[INFO] threshold={threshold}, filtered_rows={len(filtered_df)}, output_ids={len(unique_ids)}")


if __name__ == "__main__":
    start_time = time.time()
    args = parse_args()
    tmp_result_path = f"tmp/{args.tmp_name}/"
    if not os.path.exists(tmp_result_path):
        os.makedirs(tmp_result_path)

    # 读取相关性数据和种子节点数据
    correlation_df = pd.read_csv(
        args.correlation_file,
        dtype={
            "Substance 1": "string",
            "Substance 2": "string",
        },
    )

    seednode_df = pd.read_csv(
        args.seednode_file,
        dtype={"ID": "string"},
    )

    # 过滤高相关性化合物
    filter_filename = (
        tmp_result_path + "corr_pval_final_miniCor_p0.05_connectwith_seednode.csv"
    )
    filtered_df = filter_according_to_substance(
        correlation_df, seednode_df, filter_filename
    )

    # 生成高相关性节点
    high_correlated_nodes_filename = tmp_result_path + "high_correlated_nodes.csv"
    generate_high_correlated_nodes(
        filtered_df, high_correlated_nodes_filename, args.correlation_threshold
    )

    print("Spend time: ", time.time() - start_time)
