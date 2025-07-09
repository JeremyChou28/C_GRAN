import os
import argparse
import pandas as pd
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(description="Search candidates.")
    parser.add_argument(
        "--seednode_file",
        default="seednode.csv",
        type=str,
        required=True,
        help="the original seednode file containing ID and SMILES",
    )
    parser.add_argument(
        "--round_num",
        type=int,
        default=0,
        help="the round of annotation, used to name the output file",
    )
    return parser.parse_args()


def pick_final_annotation(cfmid_prediction_folder, output_file):
    # 初始化结果列表
    result_rows = []

    # 获取所有CSV文件
    csv_files = [f for f in os.listdir(cfmid_prediction_folder) if f.endswith(".csv")]

    # 遍历每个文件
    with tqdm(total=len(csv_files), desc="Processing CSV files") as pbar:
        for file in csv_files:
            file_path = os.path.join(cfmid_prediction_folder, file)
            try:
                df = pd.read_csv(file_path)

                # 跳过无CFM-ID_score列的文件
                if "CFM-ID_score" not in df.columns:
                    tqdm.write(f"Skipped {file}: No 'CFM-ID_score' column.")
                    pbar.update(1)
                    continue

                # 找出CFM-ID_score最大值的所有行
                max_cfm = df["CFM-ID_score"].max()
                candidates = df[df["CFM-ID_score"] == max_cfm]

                # 判断有无weighted_score列
                if "weighted_score" in candidates.columns:
                    # 优先按 weighted_score 选最大
                    max_weighted = candidates["weighted_score"].max()
                    candidates = candidates[
                        candidates["weighted_score"] == max_weighted
                    ]
                elif "score" in candidates.columns:
                    # 否则按 score 选最大
                    max_score = candidates["score"].max()
                    candidates = candidates[candidates["score"] == max_score]

                # 再次并列，默认选择第一行（靠前）
                best_row = candidates.iloc[0].copy()
                best_row["ID"] = os.path.splitext(file)[0]  # ✅ 去除 .csv 后缀

                # 添加到结果列表
                result_rows.append(best_row)

            except Exception as e:
                tqdm.write(f"Error processing {file}: {str(e)}")
            pbar.update(1)

    # 合并并输出
    if result_rows:
        final_df = pd.DataFrame(result_rows)
        # 将 ID 列置于最前
        cols = ["ID"] + [col for col in final_df.columns if col != "ID"]
        final_df = final_df[cols]
        final_df.to_csv(output_file, index=False)
        return final_df
    else:
        return None


def filter_annotation(df, output_file):
    # 初始化进度条
    with tqdm(total=2, desc="Filtering and selecting columns") as pbar:
        df.rename(columns={"Seednode": "Seed Node"}, inplace=True)
        df.rename(columns={"CFM-ID_score": "Score"}, inplace=True)
        df.rename(columns={"MW": "MonoIsotopic Weight"}, inplace=True)

        # Step 2: 只保留指定列
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
        df_final["ID"] = df_final["ID"].astype(int)
        # 将CID列改成超链接
        df_final["CID"] = "https://pubchem.ncbi.nlm.nih.gov/compound/" + df_final[
            "CID"
        ].astype(str)

        df_final["Seed Node"] = df_final["Seed Node"].apply(
            lambda x: "" if pd.isna(x) else str(int(x))
        )

        pbar.update(1)

    # 保存结果
    df_final.sort_values(by="Score", inplace=True)
    df_final.to_csv(output_file, index=False)
    return df_final


if __name__ == "__main__":
    args = parse_args()

    tmp_result_path = "tmp/"
    if not os.path.exists(tmp_result_path):
        os.makedirs(tmp_result_path)

    cfmid_prediction_folder = tmp_result_path + "cfmid_score_results"
    annotation_result_file = tmp_result_path + "annotation_results.csv"
    annotation_result_df = pick_final_annotation(
        cfmid_prediction_folder, annotation_result_file
    )
    if annotation_result_df is None:
        annotated_nodes_df = pd.DataFrame(columns=["ID", "SMILES"])
    else:
        output_file = tmp_result_path + "annotation_results_filtered.csv"
        annotated_nodes_df_filtered = filter_annotation(
            annotation_result_df, output_file
        )
        annotated_nodes_df = annotated_nodes_df_filtered[["ID", "SMILES"]]

    # 生成下一轮的seednode file
    seednode_df = pd.read_csv(args.seednode_file)[["ID", "SMILES"]]
    merged_df = pd.concat([seednode_df, annotated_nodes_df], ignore_index=True)
    merged_df.to_csv(
        tmp_result_path + f"annotation_with_cfmid_seednode_round{args.round_num}.csv",
        index=False,
    )
