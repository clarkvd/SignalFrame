import os
import pyBigWig
import pandas as pd
import numpy as np
from typing import List
from tqdm import tqdm


def merge_intervals_with_indices(intervals, slack):
        if not intervals:
            return []
        merged = []
        cur_start, cur_end, cur_idxs = intervals[0][0], intervals[0][1], [intervals[0][2]]
        for start, end, idx in intervals[1:]:
            if start <= cur_end + slack:
                cur_end = max(cur_end, end)
                cur_idxs.append(idx)
            else:
                merged.append((cur_start, cur_end, cur_idxs))
                cur_start, cur_end, cur_idxs = start, end, [idx]
        merged.append((cur_start, cur_end, cur_idxs))
        return merged

def smart_merge_slack(starts: np.ndarray, ends: np.ndarray) -> int:
    if len(starts) <= 1:
        return 0

    gaps = starts[1:] - ends[:-1]
    gaps = gaps[gaps > 0]
    if len(gaps) == 0:
        return 50

    # Use faster percentile approximation (interpolation="midpoint" is default, fast)
    q50 = np.percentile(gaps, 50)
    return int(min(max(q50 * 0.5, 10), 500))


def compute_signal(bigwig_path: str, bed_df: pd.DataFrame, method: str) -> pd.DataFrame:
    """
    Compute signal values from a BigWig file for each region in a BED DataFrame.

    Supports various summary statistics across each interval, such as AUC, mean, max, etc.
    Automatically merges nearby regions to reduce redundant I/O and speeds up access.

    Parameters:
    - bigwig_path (str): Path to the BigWig file.
    - bed_df (pd.DataFrame): DataFrame with BED-format intervals. Must contain 'chr', 'start', and 'end' columns.
    - method (str): Signal computation method. Must be one of:
        - 'auc': Area under the curve (sum of signal values)
        - 'mean': Average signal value
        - 'max': Maximum signal value
        - 'min': Minimum signal value
        - 'median': Median signal value
        - 'std': Standard deviation of signal values
        - 'coverage': Number of positions with non-zero signal
        - 'nonzero_mean': Mean of non-zero signal values

    Returns:
    - pd.DataFrame: Input DataFrame with an added column named after the selected `method`.

    Raises:
    - ValueError: If an invalid method is provided or if the DataFrame is missing required columns.
    """
    valid_methods = {"auc", "mean", "max", "min", "median", "std", "coverage", "nonzero_mean"}
    if method not in valid_methods:
        raise ValueError(f"Invalid method '{method}'. Choose from {valid_methods}")

    required_cols = {"chr", "start", "end"}
    if not required_cols.issubset(bed_df.columns):
        raise ValueError(f"Input DataFrame must contain columns: {required_cols}")

    bed_df = bed_df.sort_values(["chr", "start"]).copy()
    results = np.zeros(len(bed_df), dtype=np.float32)


    bw = pyBigWig.open(bigwig_path)
    try:
        chrom_sizes = bw.chroms()
        pbar = tqdm(total=len(bed_df), desc=f"Computing {method.upper()}")

        for chrom, group in bed_df.groupby("chr"):
            if chrom not in chrom_sizes:
                pbar.update(len(group))
                continue

            max_end = chrom_sizes[chrom]
            starts = group["start"].to_numpy()
            ends = group["end"].to_numpy()
            indices = group.index.to_numpy()

            slack = smart_merge_slack(starts, ends)
            intervals = list(zip(starts, ends, indices))
            merged_regions = merge_intervals_with_indices(intervals, slack=slack)

            for mstart, mend, idxs in merged_regions:
                s_clip = max(0, mstart)
                e_clip = min(mend, max_end)
                if s_clip >= e_clip:
                    for idx in idxs:
                        results[idx] = 0
                        pbar.update(1)
                    continue

                try:
                    signal = bw.values(chrom, s_clip, e_clip, numpy=True)
                except RuntimeError:
                    for idx in idxs:
                        results[idx] = 0
                        pbar.update(1)
                    continue

                signal = np.nan_to_num(signal, nan=0.0)

                for idx in idxs:
                    rstart, rend = bed_df.at[idx, "start"], bed_df.at[idx, "end"]
                    sub_start = max(rstart, s_clip)
                    sub_end = min(rend, e_clip)
                    if sub_start >= sub_end:
                        results[idx] = 0
                        pbar.update(1)
                        continue
                    offset = sub_start - s_clip
                    region = signal[offset:offset + (sub_end - sub_start)]

                    results[idx] = (
                        np.sum(region) if method == "auc" else
                        np.mean(region) if method == "mean" else
                        np.max(region) if method == "max" else
                        np.min(region) if method == "min" else
                        np.median(region) if method == "median" else
                        np.std(region) if method == "std" else
                        np.count_nonzero(region) if method == "coverage" else
                        np.mean(region[region > 0]) if method == "nonzero_mean" and np.any(region > 0) else 0
                    )
                    pbar.update(1)
        pbar.close()
    finally:
        bw.close()

    bed_df[method] = results
    return bed_df

def compute_signal_multi(bigwig_paths: List[str],
                         bed_df: pd.DataFrame,
                         method: str) -> pd.DataFrame:
    """
    Compute signal values from multiple BigWig files for each region in a BED DataFrame.

    For each input BigWig file, the function calculates a summary statistic (e.g., AUC, mean)
    over each genomic interval provided in the BED-format DataFrame. Results are added as
    new columns, one per BigWig file.

    Parameters:
    - bigwig_paths (List[str]): List of paths to BigWig files.
    - bed_df (pd.DataFrame): DataFrame with genomic intervals; must contain 'chr', 'start', and 'end' columns.
    - method (str): Signal computation method. Must be one of:
        - 'auc': Area under the curve (sum of signal values)
        - 'mean': Average signal value
        - 'max': Maximum signal value
        - 'min': Minimum signal value
        - 'median': Median signal value
        - 'std': Standard deviation of signal values
        - 'coverage': Number of positions with non-zero signal
        - 'nonzero_mean': Mean of non-zero signal values

    Returns:
    - pd.DataFrame: Input DataFrame with additional columns, one per BigWig file and method combination.

    Raises:
    - ValueError: If the input method is invalid or if required BED columns are missing.
    - ValueError: If a file is not a valid BigWig.
    """
    valid_methods = {"auc", "mean", "max", "min", "median", "std", "coverage", "nonzero_mean"}
    if method not in valid_methods:
        raise ValueError(f"Invalid method '{method}'. Choose from {valid_methods}")

    required_cols = {"chr", "start", "end"}
    if not required_cols.issubset(bed_df.columns):
        raise ValueError(f"Input DataFrame must contain columns: {required_cols}")

    bed_df = bed_df.sort_values(["chr", "start"]).copy()
    signal_matrix = np.zeros((len(bed_df), len(bigwig_paths)), dtype=np.float32)

    for bw_idx, bw_path in enumerate(bigwig_paths):
        bw = pyBigWig.open(bw_path)
        if not bw.isBigWig():
            raise ValueError(f"{bw_path} is not a valid BigWig file")

        chrom_sizes = bw.chroms()
        pbar = tqdm(total=len(bed_df), desc=f"{method.upper()}: {os.path.basename(bw_path)}")

        for chrom, group in bed_df.groupby("chr"):
            if chrom not in chrom_sizes:
                pbar.update(len(group))
                continue

            max_end = chrom_sizes[chrom]
            starts = group["start"].to_numpy()
            ends = group["end"].to_numpy()
            indices = group.index.to_numpy()

            slack = smart_merge_slack(starts, ends)
            intervals = list(zip(starts, ends, indices))
            merged_regions = merge_intervals_with_indices(intervals, slack=slack)

            for mstart, mend, idxs in merged_regions:
                s_clip = max(0, mstart)
                e_clip = min(mend, max_end)

                if s_clip >= e_clip:
                    for idx in idxs:
                        signal_matrix[idx, bw_idx] = 0
                        pbar.update(1)
                    continue

                try:
                    signal = bw.values(chrom, s_clip, e_clip, numpy=True)
                except RuntimeError:
                    for idx in idxs:
                        signal_matrix[idx, bw_idx] = 0
                        pbar.update(1)
                    continue

                signal = np.nan_to_num(signal, nan=0.0)

                for idx in idxs:
                    rstart, rend = bed_df.at[idx, "start"], bed_df.at[idx, "end"]
                    sub_start = max(rstart, s_clip)
                    sub_end = min(rend, e_clip)
                    if sub_start >= sub_end:
                        signal_matrix[idx, bw_idx] = 0
                        pbar.update(1)
                        continue

                    offset = sub_start - s_clip
                    region = signal[offset:offset + (sub_end - sub_start)]

                    signal_matrix[idx, bw_idx] = (
                        np.sum(region) if method == "auc" else
                        np.mean(region) if method == "mean" else
                        np.max(region) if method == "max" else
                        np.min(region) if method == "min" else
                        np.median(region) if method == "median" else
                        np.std(region) if method == "std" else
                        np.count_nonzero(region) if method == "coverage" else
                        np.mean(region[region > 0]) if method == "nonzero_mean" and np.any(region > 0) else 0
                    )
                    pbar.update(1)

        pbar.close()
        bw.close()

    signal_df = pd.DataFrame(
        signal_matrix,
        columns=[
            f"{os.path.splitext(os.path.basename(p))[0]}_{method}"
            for p in bigwig_paths
        ]
    )

    return pd.concat([bed_df.reset_index(drop=True), signal_df], axis=1)