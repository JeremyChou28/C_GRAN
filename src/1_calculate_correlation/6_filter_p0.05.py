import pandas as pd
from tqdm import tqdm

# 读取CSV文件
df = pd.read_csv('corr_pval_final_CD_sediment_pos_3SD_20240828_miniCor.csv')

# 筛选条件：P-Value小于0.05 
filtered_rows = []
total_rows = len(df)
with tqdm(total=total_rows) as pbar:
    for index, row in df.iterrows():
        if (row['P-Value'] < 0.05):
            filtered_rows.append(row)
        pbar.update(1)

filtered_df = pd.DataFrame(filtered_rows)

# 将筛选后的结果保存到新的CSV文件中
filtered_df.to_csv('corr_pval_final_CD_sediment_pos_3SD_20240828_miniCor_p0.05.csv', index=False)


