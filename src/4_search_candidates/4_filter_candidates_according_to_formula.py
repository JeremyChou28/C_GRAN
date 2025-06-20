import os
import pandas as pd
import re
from tqdm import tqdm

# 输入文件夹路径和输出文件夹路径
input_folder = 'output/search_pubchem_mw_result_row_level'     # 替换为你的原始CSV文件夹路径
output_folder = 'output/search_pubchem_mw_result_row_level_filtered_formula'

# 创建输出目录（如不存在）
os.makedirs(output_folder, exist_ok=True)

# 允许的元素集合（注意顺序处理避免误识别Cl为C和l）
allowed_elements = {'C', 'H', 'O', 'N', 'P', 'S', 'F', 'Cl', 'Br', 'I'}

# 用于正则提取分子式中的元素符号（如 C, Cl, Br）
element_pattern = re.compile(r'Cl|Br|[A-Z][a-z]?')

# 获取所有csv文件名
csv_files = [f for f in os.listdir(input_folder) if f.endswith('.csv')]

# 初始化进度条
with tqdm(total=len(csv_files), desc="Filtering formulas") as pbar:
    for file in csv_files:
        # 读取文件
        file_path = os.path.join(input_folder, file)
        df = pd.read_csv(file_path)

        # 检查是否有Formula列
        if 'Formula' not in df.columns:
            tqdm.write(f"Skipped {file}: No 'Formula' column found.")
            pbar.update(1)
            continue

        # 过滤非法分子式
        def is_valid_formula(formula):
            elements = element_pattern.findall(str(formula))
            return all(e in allowed_elements for e in elements)

        df_filtered = df[df['Formula'].apply(is_valid_formula)]

        # 写入输出文件
        output_path = os.path.join(output_folder, file)
        df_filtered.to_csv(output_path, index=False)

        pbar.update(1)
