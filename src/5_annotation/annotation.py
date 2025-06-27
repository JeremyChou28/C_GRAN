import os
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
        "--pubchem_mw_folder",
        type=str,
        required=True,
        help="pubchem molecular weight folder containing TTL files",
    )
    parser.add_argument(
        "--pubchem_smiles_folder",
        type=str,
        required=True,
        help="pubchem smiles folder containing TTL files",
    )
    parser.add_argument(
        "--pubchem_formula_folder",
        type=str,
        required=True,
        help="pubchem formula folder containing TTL files",
    )
    parser.add_argument(
        "--candidates_folder",
        type=str,
        help="the candidates folder to save the search results",
    )
    return parser.parse_args()



def generate_seednode_and_targetnode(molecular_network_df, seednode_df, output_filename):
    # 提取种子物质 ID 集合
    seed_ids = set(seednode_df['ID'].astype(int))  # 确保为整数类型以匹配

    # 初始化列表用于存储结果
    results = []

    # 使用 tqdm 包装 DataFrame 行迭代器
    for _, row in tqdm(molecular_network_df.iterrows(), total=len(molecular_network_df), desc="Filtering connected nodes"):
        source = int(row['source'])
        target = int(row['target'])
        corr = row['Correlation']  # 提取 Correlation 值

        source_in_seed = source in seed_ids
        target_in_seed = target in seed_ids

        # 仅保留一端在种子集合中，另一端不在
        if source_in_seed and not target_in_seed:
            results.append({
                'Seednode': source,
                'Targetnode': f"{target:.1f}",
                'Correlation': corr
            })
        elif target_in_seed and not source_in_seed:
            results.append({
                'Seednode': target,
                'Targetnode': f"{source:.1f}",
                'Correlation': corr
            })

    # 转换为 DataFrame 并按 Seednode 升序排序
    result_df = pd.DataFrame(results)
    result_df.sort_values(by='Seednode', inplace=True)

    # 保存输出文件
    result_df.to_csv(output_filename, index=False)
    return result_df


def calculate_tanimoto_similarity(seednode_df, edges_df, candidates_folder, output_folder):
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
    for _, row in tqdm(edges_df.iterrows(), total=len(edges_df), desc="Processing Edges"):
        target_id, seed_id = row["Targetnode"], row["Seednode"]
        
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
        df_target["score"] = df_target["SMILES"].apply(lambda x: tanimoto_similarity(x, seed_smiles))

        # 按 score 降序排序
        df_target = df_target.sort_values(by="score", ascending=False)
        
        # 保存到新文件
        output_filename = os.path.join(output_folder, f"{target_id}_{seed_id}.csv")
        df_target.to_csv(output_filename, index=False)

    print("相似性计算完成，结果已保存到 output 目录。")


def check_seednode(seed_target_df, output_unique, output_repeated):
    # 使用 tqdm 展示拆分进度
    with tqdm(total=1, desc="Splitting based on Targetnode count") as pbar:
        # 计算每个 Targetnode 出现次数
        target_counts = seed_target_df['Targetnode'].value_counts()

        # 获取只出现一次和多次的 Targetnode ID 集合
        unique_targets = set(target_counts[target_counts == 1].index)
        repeated_targets = set(target_counts[target_counts > 1].index)

        # 拆分为两个 DataFrame
        df_unique = seed_target_df[seed_target_df['Targetnode'].isin(unique_targets)].copy()
        df_repeated = seed_target_df[seed_target_df['Targetnode'].isin(repeated_targets)].copy()

        pbar.update(1)

    # 保存两个文件
    df_unique.to_csv(output_unique, index=False)
    df_repeated.to_csv(output_repeated, index=False)


def split_seednode():
    # 文件路径配置（请修改为你自己的路径）
    unique_csv = 'output/Seednode_and_Targetnode_unique_seednode.csv'
    not_unique_csv = 'output/Seednode_and_Targetnode_not_unique_seednodes.csv'
    source_folder = 'output/Seednode_and_Targetnode_Morgan_Similarity_score'  # 包含所有8579.0_27.0格式文件的文件夹
    output_folder_unique = 'output/Seednode_and_Targetnode_Morgan_Similarity_score_split_unique'
    output_folder_not_unique = 'output/Seednode_and_Targetnode_Morgan_Similarity_score_split_not_unique'

    # 创建输出文件夹
    os.makedirs(output_folder_unique, exist_ok=True)
    os.makedirs(output_folder_not_unique, exist_ok=True)

    # 读取两个csv文件
    df_unique = pd.read_csv(unique_csv)
    df_not_unique = pd.read_csv(not_unique_csv)

    # 构建命名集合（文件名：'Targetnode_Seednode.csv'）
    def build_filename_set(df):
        return set(f"{row['Targetnode']}_{row['Seednode']}.csv" for _, row in df.iterrows())

    unique_files = build_filename_set(df_unique)
    not_unique_files = build_filename_set(df_not_unique)

    # 遍历原始文件夹中的所有csv文件，按匹配规则分类拷贝
    all_files = [f for f in os.listdir(source_folder) if f.endswith('.csv')]

    with tqdm(total=len(all_files), desc="Classifying CSV files") as pbar:
        for file in all_files:
            source_path = os.path.join(source_folder, file)

            if file in unique_files:
                shutil.copy(source_path, os.path.join(output_folder_unique, file))
            elif file in not_unique_files:
                shutil.copy(source_path, os.path.join(output_folder_not_unique, file))
            # 忽略未匹配文件
            pbar.update(1)


if __name__ == "__main__":
    args=parse_args()
    
    tmp_result_path='tmp/'
    if not os.path.exists(tmp_result_path):
        os.makedirs(tmp_result_path)
    
    
    # 获取seed node和target node
    molecular_network_df = pd.read_csv(args.molecular_network_file)
    seednode_df = pd.read_csv(args.seednode_file)  
    seed_target_node_df=generate_seednode_and_targetnode(molecular_network_df, seednode_df, os.path.join(tmp_result_path, "Seednode_and_Targetnode.csv"))
    
    # 计算tanimoto similarity
    edges_df=pd.read_csv(args.edges_file)
    calculate_tanimoto_similarity(seednode_df, edges_df, args.candidates_folder, os.path.join(tmp_result_path, "Seednode_and_Targetnode_Morgan_Similarity_score"))
    
    # 检查种子节点
    seed_target_df=pd.read_csv(os.path.join(tmp_result_path, "Seednode_and_Targetnode.csv"))
    output_unique = tmp_result_path+'Seednode_and_Targetnode_unique_seednode.csv'
    output_repeated = tmp_result_path+'Seednode_and_Targetnode_not_unique_seednodes.csv'
    check_seednode(seed_target_df, output_unique, output_repeated)