import numpy as np
import os
import pandas as pd
import argparse
from tqdm import tqdm
import multiprocessing
from multiprocessing import shared_memory, Process, Pool
from functools import partial
from scipy.stats import pearsonr
import pickle

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_data",
        help="input data path, TXT file",
        type=str
    )
    parser.add_argument(
        "--filter-expt",
        default=False,
        help="filter out exception values",
        action="store_true",
    )
    parser.add_argument(
        "--need-transpose",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--n-jobs",
        default=1,
        type=int,
        help="number of jobs",
    )
    parser.add_argument(
        "--run-name",
        default="",
        type=str,
        help="run name for output pickle filename",
    )
    return parser.parse_args()

def filter_out_exceptions(row: pd.Series):
    q1 = row.quantile(0.25)
    q3 = row.quantile(0.75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    row[(row < lower_bound) | (row > upper_bound)] = np.nan
    return row

def pearson_worker(row_j, row_i):
    mask_i = ~np.isnan(row_i)
    mask_j = ~np.isnan(row_j)
    mask_valid = mask_i & mask_j
    common_obs = np.sum(mask_valid)
    if common_obs < 3:
        return (np.nan, np.nan, 0)
    corr, p_value = pearsonr(row_i[mask_valid], row_j[mask_valid])
    return corr, p_value, common_obs

if __name__ == "__main__":
    args = parse_args()

    print(f'Loading data {args.input_data}...')
    df = pd.read_csv(args.input_data, sep='\t')  
    names = df.columns
    if args.need_transpose:
        df = df.T
    assert df.shape[1] == 60, "the number of columns should be 60, add or delete the --need-transpose flag"  #需要由用户输入样品数量
    print(f'Data loaded. shape: {df.shape}')

    if args.filter_expt:
        print("filter expt values")
        for row_id, row in df.iterrows():
            filter_out_exceptions(row)
    else:
        print("not filter expt values")
    out_file_name = f"corr_pval_{args.run_name}.pickle"
    print(f"out_file_name: {out_file_name}")

    arr = df.to_numpy()
    tot_iters = len(df) * (len(df) - 1) // 2
    pbar = tqdm(total=tot_iters)
    res_shape = (arr.shape[0], arr.shape[0])
    corr_arr = np.zeros(res_shape)
    pval_arr = np.zeros(res_shape)
    non_zero_count_arr = np.zeros(res_shape, dtype=int)

    print("number of jobs:", args.n_jobs)
    with multiprocessing.Pool(processes=args.n_jobs) as pool:
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
    
    with open(out_file_name, "wb") as f:
        pickle.dump(save_res, f)
