import pandas as pd
import numpy as np
from typing import List, Union, Optional

def normalize_signal(df: pd.DataFrame,
                     columns: Union[str, List[str]],
                     method: str,
                     pseudocount: float = 0.1,
                     reference_matrix: Optional[Union[str, np.ndarray]] = None) -> pd.DataFrame:
    """
    Normalize one or more signal columns using common bioinformatics transformations.

    Parameters:
    - df (pd.DataFrame): Input DataFrame with signal values.
    - columns (str or List[str]): One or more column names to normalize.
    - method (str): Normalization method. Options:
        - 'length': Normalize by region length (AUC per bp).
        - 'zscore': Z-score normalization.
        - 'log2': Log2(x + pseudocount).
        - 'minmax': Scale values to [0, 1].
        - 'quantile': Match rank distribution to a reference.
            - If reference_matrix == 'median', use the median profile.
            - If reference_matrix == 'mean', use the mean profile.
    - pseudocount (float): Value added before log2 transformation (default = 0.1).
    - reference_matrix (str or np.ndarray): Used only for 'quantile'. If str, must be 'median' or 'mean'.
        If ndarray, must have shape (n_rows, n_columns).

    Returns:
    - pd.DataFrame: New DataFrame with normalized columns appended as '<column>_norm'.

    Raises:
    - ValueError: On invalid column names, methods, or shape mismatches.
    """
    df = df.copy()
    if isinstance(columns, str):
        columns = [columns]

    # Compute reference matrix if needed
    if method == "quantile" and isinstance(reference_matrix, str):
        values = df[columns].values
        sorted_values = np.sort(values, axis=0)
        if reference_matrix.lower() == "median":
            ref = np.median(sorted_values, axis=1).reshape(-1, 1)
        elif reference_matrix.lower() == "mean":
            ref = np.mean(sorted_values, axis=1).reshape(-1, 1)
        else:
            raise ValueError("reference_matrix must be 'median', 'mean', or an ndarray")
        reference_matrix = np.repeat(ref, len(columns), axis=1)

    for i, col in enumerate(columns):
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in DataFrame")

        if method == "length":
            if "start" not in df.columns or "end" not in df.columns:
                raise ValueError("Columns 'start' and 'end' required for length normalization")
            df[col + "_norm"] = df[col] / (df["end"] - df["start"])

        elif method == "zscore":
            mean = df[col].mean()
            std = df[col].std()
            df[col + "_norm"] = (df[col] - mean) / std if std > 0 else 0.0

        elif method == "log2":
            df[col + "_norm"] = np.log2(df[col] + pseudocount)

        elif method == "minmax":
            min_val = df[col].min()
            max_val = df[col].max()
            df[col + "_norm"] = 0.0 if max_val == min_val else (df[col] - min_val) / (max_val - min_val)

        elif method == "quantile":
            if reference_matrix is None:
                raise ValueError("reference_matrix is required for quantile normalization")
            if reference_matrix.shape != (len(df), len(columns)):
                raise ValueError(f"reference_matrix must have shape ({len(df)}, {len(columns)})")

            ref_col = reference_matrix[:, i]
            ranks = np.argsort(np.argsort(df[col].values))
            ref_sorted = np.sort(ref_col)
            df[col + "_norm"] = ref_sorted[ranks]

        else:
            raise ValueError("method must be one of: 'length', 'zscore', 'log2', 'minmax', 'quantile'")

    return df

def compare_tracks(df: pd.DataFrame,
                   reference: str,
                   comparisons: Union[str, List[str]],
                   mode: Union[str, List[str]] = ["difference", "fold_change", "log2FC", "percent_change"],
                   pseudocount: float = 0.1) -> pd.DataFrame:
    """
    Compare signal tracks per region between a reference and one or more other tracks.

    Parameters:
    - df (pd.DataFrame): DataFrame with signal values.
    - reference (str): Column name to use as the reference.
    - comparisons (str or list of str): Columns to compare against the reference.
    - mode (str or list of str): Comparison types. Options:
        - 'difference': ref - target
        - 'fold_change': (ref + pseudocount) / (target + pseudocount)
        - 'log2FC': log2((ref + pseudocount) / (target + pseudocount))
        - 'percent_change': (ref - target) / (target + pseudocount)
    - pseudocount (float): Value to prevent divide-by-zero and log issues.

    Returns:
    - pd.DataFrame: DataFrame with added comparison columns.
    """
    if isinstance(comparisons, str):
        comparisons = [comparisons]
    if isinstance(mode, str):
        mode = [mode]

    if reference not in df.columns:
        raise ValueError(f"Reference column '{reference}' not found in DataFrame")
    for target in comparisons:
        if target not in df.columns:
            raise ValueError(f"Comparison column '{target}' not found in DataFrame")

    df = df.copy()

    for target in comparisons:
        for mode in mode:
            if mode == "difference":
                df[f"{reference}_vs_{target}_diff"] = df[reference] - df[target]

            elif mode == "fold_change":
                df[f"{reference}_vs_{target}_FC"] = (
                    (df[reference] + pseudocount) / (df[target] + pseudocount)
                )

            elif mode == "log2FC":
                df[f"{reference}_vs_{target}_log2FC"] = np.log2(
                    (df[reference] + pseudocount) / (df[target] + pseudocount)
                )

            elif mode == "percent_change":
                df[f"{reference}_vs_{target}_pct_change"] = (
                    (df[reference] - df[target]) / (df[target] + pseudocount)
                )

            else:
                raise ValueError(f"Invalid comparison mode: {mode}")

    return df

def sort_signal_df(df: pd.DataFrame,
                sort_by: str = "genomic_position",
                ascending: bool = True) -> pd.DataFrame:
    """
    Sort AUC result DataFrame by genomic position or a specified column.

    Parameters:
    - df (pd.DataFrame): DataFrame with AUC and genomic information.
    - sort_by (str): Sorting method:
        - 'genomic_position': Sort by 'chr' and 'start'.
        - Any column name in the DataFrame.
    - ascending (bool): Sort order. True for ascending (default), False for descending.

    Returns:
    - pd.DataFrame: Sorted copy of the original DataFrame.

    Raises:
    - ValueError: If 'sort_by' is not 'genomic_position' or a valid column name.
    """
    if sort_by == "genomic_position":
        if not {"chr", "start"}.issubset(df.columns):
            raise ValueError("DataFrame must contain 'chr' and 'start' columns for genomic sorting.")
        return df.sort_values(by=["chr", "start"], ascending=ascending).reset_index(drop=True)

    elif sort_by in df.columns:
        return df.sort_values(by=sort_by, ascending=ascending).reset_index(drop=True)

    else:
        raise ValueError(
            f"Invalid sort_by value: '{sort_by}'. Must be 'genomic_position' or a column in the DataFrame."
        )