import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, DataStructs
from tqdm import tqdm
import warnings
# 忽略警告信息
warnings.filterwarnings("ignore")

# 文件路径
seednode_file = "test_files/14seednode.csv"
edges_file = "tmp/Seednode_and_Targetnode.csv"
targetnode_dir = "../4_search_candidates/candidates"  # 目标物质的 CSV 文件存放目录
output_dir = "tmp/Seednode_and_Targetnode_Morgan_Similarity_score"  # 结果存放目录
os.makedirs(output_dir, exist_ok=True)

# 读取 Seednode 的 SMILES 信息
df_seednode = pd.read_csv(seednode_file)
seednode_smiles_dict = dict(zip(df_seednode["ID"], df_seednode["SMILES"]))

# 读取 Targetnode 和 Seednode 关系对
df_edges = pd.read_csv(edges_file)

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
for _, row in tqdm(df_edges.iterrows(), total=len(df_edges), desc="Processing Edges"):
    target_id, seed_id = int(row["Targetnode"]), int(row["Seednode"])
    
    # 获取 Seednode 的 SMILES
    seed_smiles = seednode_smiles_dict.get(seed_id)
    if not seed_smiles:
        continue  # 如果 Seednode 没有 SMILES，跳过
    
    # 读取 Targetnode 文件
    target_file = os.path.join(targetnode_dir, f"{target_id}.csv")
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
    output_filename = os.path.join(output_dir, f"{target_id}_{seed_id}.csv")
    df_target.to_csv(output_filename, index=False)

print("相似性计算完成，结果已保存到 output 目录。")
