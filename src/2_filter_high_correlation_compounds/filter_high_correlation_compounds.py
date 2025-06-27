import os
import time
import pandas as pd
from tqdm import tqdm
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Filter high correlation compounds and generate high correlated nodes.")
    parser.add_argument(
        "--seednode_file",
        default="14seednode.csv",
        type=str,
        required=True,
        help="input seednode file",
    )
    parser.add_argument(
        "--correlation_file",
        type=str,
        required=True,
        help="input correlation file",
    )
    return parser.parse_args()

def filter_according_to_substance(correlation_df, seednode_df, output_filename):
    # 提取第二个CSV文件中的ID列
    substances = seednode_df['ID'].tolist()

    # 筛选Substance 1或者Substance 2在ID中的行
    filtered_df = correlation_df[(correlation_df['Substance 1'].isin(substances)) | (correlation_df['Substance 2'].isin(substances))]

    # 保存筛选后的结果到新的CSV文件，添加进度条
    with tqdm(total=filtered_df.shape[0]) as pbar:
        filtered_df.to_csv(output_filename, index=False)
        pbar.update(filtered_df.shape[0])
    return filtered_df

def generate_high_correlated_nodes(df, output_filename):
    # 筛选出 Correlation 大于 0.7 的行  (该阈值由用户自定义)
    filtered_df = df[df['Correlation'] > 0.7]

    # 初始化进度条：总共处理两列数据
    with tqdm(total=2, desc="Processing IDs") as pbar:
        # 提取Substance 1和Substance 2列
        substances = pd.concat([filtered_df['Substance 1'], filtered_df['Substance 2']])
        pbar.update(1)

        # 去重并按升序排序
        unique_ids = pd.DataFrame(sorted(substances.unique()), columns=['ID'])
        pbar.update(1)

    # 保存为新的CSV文件
    unique_ids.to_csv(output_filename, index=False)


if __name__ == "__main__":
    start_time = time.time()
    args= parse_args()
    tmp_result_path = 'tmp/'
    if not os.path.exists(tmp_result_path):
        os.makedirs(tmp_result_path)
        
    # 读取相关性数据和种子节点数据
    correlation_df = pd.read_csv(args.correlation_file)
    seednode_df = pd.read_csv(args.seednode_file)
    
    # 过滤高相关性化合物
    filter_filename=tmp_result_path+'corr_pval_final_CD_sediment_pos_3SD_20240828_miniCor_p0.05_connectwith_14seednode.csv'
    filtered_df=filter_according_to_substance(correlation_df, seednode_df,filter_filename)
    
    # 生成高相关性节点
    high_correlated_nodes_filename=tmp_result_path+'high_correlated_nodes.csv'
    generate_high_correlated_nodes(filtered_df,high_correlated_nodes_filename)
    
    print("Spend time: ", time.time() - start_time)