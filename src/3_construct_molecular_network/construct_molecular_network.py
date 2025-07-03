import os
import time
import argparse
import pandas as pd
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(description="Construct molecular network.")
    # online parameters
    parser.add_argument(
        "--source_target_file",
        default="source_target.csv",
        type=str,
        required=True,
        help="input seednode file",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.7,
        required=True,
        help="the threshold for filtering high correlation compounds",
    )
    parser.add_argument(
        "--RT_threshold",
        type=float,
        default=0.01,
        required=True,
        help="the threshold for filtering based on retention time (RT)",
    )
    # offline parameters
    parser.add_argument(
        "--correlation_file",
        type=str,
        required=True,
        help="input correlation file",
    )
    return parser.parse_args()


def merge_source_target_cor(source_target_df, correlation_df, output_filename):
    # 转换为整数类型，确保匹配一致性
    source_target_df["source"] = source_target_df["source"].astype(int)
    source_target_df["target"] = source_target_df["target"].astype(int)
    correlation_df["Substance 1"] = correlation_df["Substance 1"].astype(int)
    correlation_df["Substance 2"] = correlation_df["Substance 2"].astype(int)

    # 构建无序配对键（使 source/target 与 Substance 1/2 顺序无关）
    source_target_df["pair_key"] = (
        source_target_df[["source", "target"]].min(axis=1).astype(str)
        + "_"
        + source_target_df[["source", "target"]].max(axis=1).astype(str)
    )
    correlation_df["pair_key"] = (
        correlation_df[["Substance 1", "Substance 2"]].min(axis=1).astype(str)
        + "_"
        + correlation_df[["Substance 1", "Substance 2"]].max(axis=1).astype(str)
    )

    # 精简 correlation_df 只保留用于合并的列
    correlation_df_slim = correlation_df[["pair_key", "Correlation", "P-Value"]]

    # 合并数据并加入进度条
    with tqdm(total=1, desc="Merging correlation data") as pbar:
        merged = pd.merge(
            source_target_df, correlation_df_slim, on="pair_key", how="left"
        )
        pbar.update(1)

    # 删除构造用的临时列
    merged.drop(columns=["pair_key"], inplace=True)

    # 保存输出文件
    merged.to_csv(output_filename, index=False)
    return merged


def edit_molecular_network(
    source_target_cor_df, output_filename, threshold, RT_threshold
):
    # 初始化进度条（2 步）
    with tqdm(total=2, desc="Filtering data") as pbar:

        # Step 1: 保留 Correlation > 0.7 的行  #这个相关性的阈值为之前筛选高相关物质时用户输入的阈值
        df_filtered = source_target_cor_df[
            source_target_cor_df["Correlation"] > threshold
        ].copy()
        pbar.update(1)

        # Step 2: 删除 abs(source_RT - target_RT) < 0.01 且 Correlation > 0.9 的行
        condition = (
            abs(df_filtered["source_RT"] - df_filtered["target_RT"]) < RT_threshold
        ) & (df_filtered["Correlation"] > 0.9)
        df_final = df_filtered[~condition].copy()
        pbar.update(1)

    # 保存结果到新文件
    df_final.to_csv(output_filename, index=False)


if __name__ == "__main__":
    start_time = time.time()
    args = parse_args()

    tmp_result_path = "tmp/"
    if not os.path.exists(tmp_result_path):
        os.makedirs(tmp_result_path)

    # 读取分子网络source node和target node
    source_target_df = pd.read_csv(args.source_target_file)

    # 读取相关性数据
    correlation_df = pd.read_csv(args.correlation_file)

    # 合并source target和相关性数据
    source_target_cor_filename = tmp_result_path + "source_target_cor.csv"
    source_target_cor_df = merge_source_target_cor(
        source_target_df, correlation_df, source_target_cor_filename
    )

    # 编辑分子网络
    output_filename = tmp_result_path + "source_target_cor_edit.csv"
    edit_molecular_network(
        source_target_cor_df, output_filename, args.threshold, args.RT_threshold
    )

    print("Spend time: ", time.time() - start_time)
