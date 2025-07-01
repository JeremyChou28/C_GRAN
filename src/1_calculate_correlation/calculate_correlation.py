import os
import time
import pickle
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import norm
import multiprocessing
from functools import partial
from scipy.stats import pearsonr
from multiprocessing import shared_memory, Process, Pool
import warnings

warnings.filterwarnings("ignore")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-data",
        default="test.txt",
        type=str,
        required=True,
        help="input data path, TXT file",
    )
    parser.add_argument(
        "--compounds-num",
        type=int,
        required=True,
        help="number of compounds")
    parser.add_argument(
        "--samples-num",
        type=int,
        required=True,
        help="number of samples")
    parser.add_argument(
        "--correlation-result-filename",
        default='correlation_results.csv',
        type=str,
        help="the calculated correlation result filename, which will be saved as a CSV file",
    )
    parser.add_argument(
        "--n-jobs",
        default=1,
        type=int,
        help="number of jobs",
    )
    parser.add_argument(
        "--run-name",
        default="test",
        type=str,
        help="run name for output pickle filename",
    )
    return parser.parse_args()

def filter_out_exceptions(df):
    for row_id, row in df.iterrows():
        q1 = row.quantile(0.25)
        q3 = row.quantile(0.75)
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        row[(row < lower_bound) | (row > upper_bound)] = np.nan
    return df

def pearson_worker(row_j, row_i):
    mask_i = ~np.isnan(row_i)
    mask_j = ~np.isnan(row_j)
    mask_valid = mask_i & mask_j
    common_obs = np.sum(mask_valid)
    if common_obs < 3:
        return (np.nan, np.nan, 0)
    corr, p_value = pearsonr(row_i[mask_valid], row_j[mask_valid])
    return corr, p_value, common_obs

def calculate_correlations(df,jobs):
    names = df.index
    
    arr = df.to_numpy()
    tot_iters = len(df) * (len(df) - 1) // 2
    pbar = tqdm(total=tot_iters)
    res_shape = (arr.shape[0], arr.shape[0])
    corr_arr = np.zeros(res_shape)
    pval_arr = np.zeros(res_shape)
    non_zero_count_arr = np.zeros(res_shape, dtype=int)

    print("number of jobs:", jobs)
    with multiprocessing.Pool(processes=jobs) as pool:
        for i, row in enumerate(arr):
            if i > 0:
                worker4i = partial(pearson_worker, row_i=arr[i])
                results = pool.map(worker4i, arr[:i])
                for j in range(i):
                    corr_arr[i][j], pval_arr[i][j], non_zero_count_arr[i][j] = results[j]
            pbar.update(i)
    
    corr_df = pd.DataFrame(corr_arr, columns=names, index=names)
    pval_df = pd.DataFrame(pval_arr, columns=names, index=names)
    non_zero_count_df = pd.DataFrame(non_zero_count_arr, columns=names, index=names)
    save_res = {
        'corr': corr_df,
        'p_value': pval_df,
        'non_zero_count': non_zero_count_df
    }
    
    return save_res


def save_to_csv(correlation_results, csv_filename):
    f_corr_df, f_pv_df , f_n_df = correlation_results['corr'], correlation_results['p_value'], correlation_results['non_zero_count']
    
    # 获取物质名称
    substances = f_corr_df.index.tolist()

    # 创建一个空的DataFrame用于存储数据
    combined_data = []

    # 计算总共需要循环的次数
    total_iterations = sum(range(len(substances)))

    # 使用tqdm来显示进度条
    with tqdm(total=total_iterations, desc="Progress") as pbar:
        # 循环遍历所有物质对，并将相关系数和P值和n值存储到新的DataFrame中
        for i, substance1 in enumerate(substances):
            for j, substance2 in enumerate(substances):
                if j < i:  # 只输出下三角部分
                    correlation = f_corr_df.iloc[i, j]
                    p_value = f_pv_df.iloc[i, j]
                    n = f_n_df.iloc[i,j]
                    combined_data.append([substance1, substance2, correlation, p_value, n])
                    if (i * len(substances) + j) % 100 == 0:  # 每100次更新一次进度
                        pbar.update(100)  # 更新进度条

    # 创建DataFrame并设置列名
    combined_df = pd.DataFrame(combined_data, columns=['Substance 1', 'Substance 2', 'Correlation', 'P-Value', 'n'])

    # 将数据保存到CSV文件中
    combined_df.to_csv(csv_filename, index=False)
    return combined_df

def compare_filtered_unfiltered(df1,df2,output_filename):
    # 检查两个DataFrame是否具有相同的结构
    if df1.shape != df2.shape:
        print("两个文件的结构不一致，无法比较。")
        exit()

    # 创建一个列表来存储不同的相关系数及P值
    different_correlations = []

    # 获取总行数并计算进度条更新的次数
    total_rows = len(df1)
    num_updates = total_rows // 100

    # 创建进度条
    progress_bar = tqdm(total=total_rows, desc="Comparing correlations")

    # 记录不同的相关系数及其位置
    counter = 0
    for index, row in df1.iterrows():
        # 检查相同位置的相关系数是否相同
        corr1 = row['Correlation']
        corr2 = df2.iloc[index]['Correlation']
        pval1 = row['P-Value']
        pval2 = df2.iloc[index]['P-Value']
        n1 = row['n']
        n2 =df2.iloc[index]['n']
        if pd.isnull(corr1) and pd.isnull(corr2):
            continue  # 如果两个相关系数都是缺失值，则跳过比较
        elif pd.isnull(corr1) or pd.isnull(corr2):
            different_correlations.append([int(row['Substance 1']), int(row['Substance 2']), corr1, pval1, n1, corr2, pval2, n2])
        elif round(corr1, 9) != round(corr2, 9):
            different_correlations.append([int(row['Substance 1']), int(row['Substance 2']), corr1, pval1, n1, corr2, pval2, n2])
        
        counter += 1
        if counter % 100 == 0:
            progress_bar.update(100)

    # 关闭进度条
    progress_bar.close()

    # 将不同的相关系数及P值写入CSV文件
    different_correlations_df = pd.DataFrame(different_correlations, columns=['Substance 1', 'Substance 2', 'Correlation 1', 'P-Value 1', 'n1', 'Correlation 2', 'P-Value 2', 'n2'])
    different_correlations_df.to_csv(output_filename, index=False)
    return different_correlations_df


def fisher_z_transform(r, n):
    # 检查相关系数是否在 (-1, 1) 区间内，避免除以零和取对数出现负数的情况
    valid_r = np.clip(r, -0.999999999, 0.999999999)
    return 0.5 * (pd.Series(r).apply(lambda x: 0.5 * (np.log(1 + x) - np.log(1 - x))) / np.sqrt(1 / (n - 3)))

def z_test(z1, z2, n1, n2):
    se = np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    z_statistic = (z1 - z2) / se
    p_value = 2 * (1 - norm.cdf(abs(z_statistic)))
    return p_value


def exec_fish_z(df,output_filename):
    # 初始化进度条
    total_iterations = len(df)
    with tqdm(total=total_iterations, desc="Processing Rows") as pbar:
        # 进行Fisher's z变换
        df['z1'] = np.where(df['n1'] <= 3, np.nan, fisher_z_transform(df['Correlation 1'], df['n1']))
        df['z2'] = np.where(df['n2'] <= 3, np.nan, fisher_z_transform(df['Correlation 2'], df['n2']))
        pbar.update(len(df))  # 更新进度条

        # 进行Z检验
        df['p_value'] = np.where((df['n1'] <= 3) | (df['n2'] <= 3), np.nan, z_test(df['z1'], df['z2'], df['n1'], df['n2']))
        pbar.update(len(df))  # 更新进度条

    # 选择显著差异的行
    significant_rows = df[(df['p_value'] < 0.05) | df['p_value'].isna()]

    # 输出到CSV文件
    significant_rows.to_csv(output_filename, index=False)
    return significant_rows

def replace_with_miniCor(unfiltered_df, significant_df, output_filename):
    
    # 使用Substance 1和Substance 2作为索引合并两个数据框
    merged_df = pd.merge(unfiltered_df, significant_df, on=['Substance 1', 'Substance 2'], how='left')

    # 判断较小的Correlation值，并获取对应的P-Value
    # 创建布尔掩码：True 表示 Correlation 1 更小，False 表示 Correlation 2 更小或相等
    mask = merged_df['Correlation 1'].abs() < merged_df['Correlation 2'].abs()

    # 创建新的列 Correlation_min 和对应的 P-Value_min
    merged_df['Correlation_min'] = merged_df['Correlation 1'].where(mask, merged_df['Correlation 2'])
    merged_df['P-Value_min'] = merged_df['P-Value 1'].where(mask, merged_df['P-Value 2'])

    # 将第一个文件中的 Correlation 和 P-Value 列替换为 Correlation_min 和 P-Value_min（若存在）
    merged_df['Correlation'] = merged_df['Correlation_min'].fillna(merged_df['Correlation'])
    merged_df['P-Value'] = merged_df['P-Value_min'].fillna(merged_df['P-Value'])

    # 删除第二个文件中不需要的列
    merged_df.drop(['Correlation 1', 'Correlation 2', 'P-Value 1', 'P-Value 2', 'Correlation_min', 'P-Value_min'], axis=1, inplace=True)

    # 保存结果到新的CSV文件，只包含需要的四列，并添加进度条
    with tqdm(total=merged_df.shape[0]) as pbar:
        merged_df[['Substance 1', 'Substance 2', 'Correlation', 'P-Value']].to_csv(output_filename, index=False)
        pbar.update(merged_df.shape[0])

    return merged_df

def filter_p(df, output_filename):
    filtered_df = df[df['P-Value'] < 0.05]
    filtered_df.to_csv(output_filename, index=False)


if __name__ == "__main__":
    start_time = time.time()
    args = parse_args()

    print(f'Loading data {args.input_data}...')
    df = pd.read_csv(args.input_data, sep='\t')  
    
    assert df.shape[0] == args.compounds_num and df.shape[1] == args.samples_num, "The shape read from the file should be (number of substances, number of samples)" 
    
    print(f'Data loaded. shape: {df.shape}')

    tmp_result_path="tmp/"
    if not os.path.exists(tmp_result_path):
        os.makedirs(tmp_result_path)

    # 过滤异常值
    filter_df=df.copy()
    filtered_expt_df=filter_out_exceptions(filter_df)
    filter_expt_correlations=calculate_correlations(filtered_expt_df, args.n_jobs)
    filtered_csv_filename = tmp_result_path+f"corr_pval_with_n_{args.run_name}_filtered.csv"
    filtered_csv_df=save_to_csv(filter_expt_correlations, filtered_csv_filename)
    
    # 保留异常值
    unfilter_df=df.copy()
    unfiltered_correlations=calculate_correlations(unfilter_df, args.n_jobs)
    unfiltered_csv_filename = tmp_result_path+f"corr_pval_with_n_{args.run_name}_unfiltered.csv"
    unfiltered_csv_df=save_to_csv(unfiltered_correlations, unfiltered_csv_filename)
    
    # 比较过滤和未过滤的结果
    compare_results_filename=tmp_result_path+f"Sediment_pos_3SD_20240812_different_correlations_with_n.csv"
    different_correlations_df=compare_filtered_unfiltered(filtered_csv_df, unfiltered_csv_df, compare_results_filename)
    
    # 进行Fisher's z变换和Z检验
    fish_z_results_filename=tmp_result_path+f"significant_Sediment_pos_3SD_20240828_different_correlations_with_n_true.csv"
    significant_rows=exec_fish_z(different_correlations_df,fish_z_results_filename)
  
    # 替换为较小的相关系数和对应的P值
    replace_with_miniCor_filename=tmp_result_path+f'corr_pval_final_CD_sediment_pos_3SD_20240828_miniCor.csv'
    replace_with_miniCor_df=replace_with_miniCor(unfiltered_csv_df, significant_rows, replace_with_miniCor_filename)
    
    # filter p值小于0.05
    filter_p(replace_with_miniCor_df, args.correlation_result_filename)
    
    print("Spend time: ", time.time() - start_time)