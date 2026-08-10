"""Small-scale reference implementation copied from the pre-refactor logic."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.stats import norm, pearsonr
from statsmodels.stats.multitest import multipletests


SUBSTANCE_COLUMNS = ["Substance 1", "Substance 2"]


def normalize_substance_id_series(series):
    normalized = series.astype("string").str.strip()
    return normalized.mask(normalized == "")


def normalize_substance_columns(df, columns=SUBSTANCE_COLUMNS):
    normalized_df = df.copy()
    for column in columns:
        normalized_df[column] = normalize_substance_id_series(
            normalized_df[column]
        )
    return normalized_df


def filter_out_exceptions(df):
    # Intentionally retained in the original row-wise pandas form.
    for _, row in df.iterrows():
        q1 = row.quantile(0.25)
        q3 = row.quantile(0.75)
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        row[(row < lower_bound) | (row > upper_bound)] = np.nan
    return df


def calculate_pair_dataframe(values: np.ndarray, substance_ids) -> pd.DataFrame:
    substance_ids = [str(value).strip() for value in substance_ids]
    records = []
    for i, substance1 in enumerate(substance_ids):
        for j in range(i):
            common = ~np.isnan(values[i]) & ~np.isnan(values[j])
            n = int(common.sum())
            if n < 3:
                corr = np.nan
                p_value = np.nan
            else:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    result = pearsonr(values[i, common], values[j, common])
                corr = float(result.statistic)
                p_value = float(result.pvalue)
            records.append(
                [
                    substance1,
                    substance_ids[j],
                    corr,
                    p_value,
                    n,
                ]
            )
    return pd.DataFrame(
        records,
        columns=[
            "Substance 1",
            "Substance 2",
            "Correlation",
            "P-Value",
            "n",
        ],
    )


def compare_filtered_unfiltered(df1, df2):
    df1 = normalize_substance_columns(df1)
    df2 = normalize_substance_columns(df2)
    different_correlations = []
    for position in range(len(df1)):
        row1 = df1.iloc[position]
        row2 = df2.iloc[position]
        corr1 = row1["Correlation"]
        corr2 = row2["Correlation"]
        pval1 = row1["P-Value"]
        pval2 = row2["P-Value"]
        n1 = row1["n"]
        n2 = row2["n"]

        both_missing = pd.isna(corr1) and pd.isna(corr2)
        one_missing = pd.isna(corr1) or pd.isna(corr2)
        values_different = (
            not one_missing
            and round(float(corr1), 9) != round(float(corr2), 9)
        )
        if both_missing:
            continue
        if one_missing or values_different:
            different_correlations.append(
                [
                    row1["Substance 1"],
                    row1["Substance 2"],
                    corr1,
                    pval1,
                    n1,
                    corr2,
                    pval2,
                    n2,
                ]
            )
    return pd.DataFrame(
        different_correlations,
        columns=[
            "Substance 1",
            "Substance 2",
            "Correlation 1",
            "P-Value 1",
            "n1",
            "Correlation 2",
            "P-Value 2",
            "n2",
        ],
    )


def fisher_z_transform(r, n):
    valid_r = np.clip(r, -0.999999999, 0.999999999)
    return 0.5 * (
        pd.Series(valid_r).apply(
            lambda x: 0.5 * (np.log(1 + x) - np.log(1 - x))
        )
        / np.sqrt(1 / (n - 3))
    )


def z_test(z1, z2, n1, n2):
    se = np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    z_statistic = (z1 - z2) / se
    return 2 * (1 - norm.cdf(abs(z_statistic)))


def exec_fish_z(df):
    df = normalize_substance_columns(df)
    df["z1"] = np.where(
        df["n1"] <= 3,
        np.nan,
        fisher_z_transform(df["Correlation 1"], df["n1"]),
    )
    df["z2"] = np.where(
        df["n2"] <= 3,
        np.nan,
        fisher_z_transform(df["Correlation 2"], df["n2"]),
    )
    df["p_value"] = np.where(
        (df["n1"] <= 3) | (df["n2"] <= 3),
        np.nan,
        z_test(df["z1"], df["z2"], df["n1"], df["n2"]),
    )
    return df[(df["p_value"] < 0.05) | df["p_value"].isna()].copy()


def replace_with_mini_cor(unfiltered_df, significant_df):
    unfiltered_df = normalize_substance_columns(unfiltered_df)
    significant_df = normalize_substance_columns(significant_df)
    merged_df = pd.merge(
        unfiltered_df,
        significant_df,
        on=SUBSTANCE_COLUMNS,
        how="left",
        validate="one_to_one",
    )
    mask = merged_df["Correlation 1"].abs() < merged_df["Correlation 2"].abs()
    merged_df["Correlation_min"] = merged_df["Correlation 1"].where(
        mask,
        merged_df["Correlation 2"],
    )
    merged_df["P-Value_min"] = merged_df["P-Value 1"].where(
        mask,
        merged_df["P-Value 2"],
    )
    merged_df["Correlation"] = merged_df["Correlation_min"].fillna(
        merged_df["Correlation"]
    )
    merged_df["P-Value"] = merged_df["P-Value_min"].fillna(
        merged_df["P-Value"]
    )
    return merged_df[[
        "Substance 1",
        "Substance 2",
        "Correlation",
        "P-Value",
    ]].copy()


def run_legacy_pipeline(values: np.ndarray, substance_ids):
    frame = pd.DataFrame(values.copy(), index=substance_ids)
    filtered = filter_out_exceptions(frame.copy()).to_numpy(dtype=np.float64)
    filtered_pairs = calculate_pair_dataframe(filtered, substance_ids)
    unfiltered_pairs = calculate_pair_dataframe(values, substance_ids)
    different = compare_filtered_unfiltered(filtered_pairs, unfiltered_pairs)
    significant_difference = exec_fish_z(different)
    before_bh = replace_with_mini_cor(
        unfiltered_pairs,
        significant_difference,
    )

    adjusted = np.full(len(before_bh), np.nan, dtype=np.float64)
    valid = before_bh["P-Value"].notna().to_numpy()
    if np.any(valid):
        adjusted[valid] = multipletests(
            before_bh.loc[valid, "P-Value"].to_numpy(),
            method="fdr_bh",
        )[1]
    with_adjusted = before_bh.copy()
    with_adjusted["adjusted p-value"] = adjusted
    final = with_adjusted[with_adjusted["adjusted p-value"] < 0.05].copy()
    return before_bh, final, filtered
