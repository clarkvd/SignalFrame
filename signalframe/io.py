import pandas as pd
from typing import List, Optional

def load_bed(bed_path: str, extra_col_names: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Load a BED file into a pandas DataFrame with standardized column names.

    Parameters:
    - bed_path (str): Path to the BED file.
    - extra_col_names (Optional[List[str]]): Optional list of names for any columns beyond the first three ('chr', 'start', 'end').
      If provided, the length must match the number of extra columns in the file.

    Returns:
    - pd.DataFrame: A DataFrame with at least the columns 'chr', 'start', and 'end'.
      Additional columns are either automatically named (col4, col5, ...) or given custom names from extra_col_names.

    Raises:
    - ValueError: If the file has fewer than 3 columns.
    - ValueError: If the number of provided extra_col_names does not match the number of extra columns.
    """
    df = pd.read_csv(bed_path, sep='\t', header=None)

    if df.shape[1] < 3:
        raise ValueError("BED file must have at least 3 columns")

    base_cols = ['chr', 'start', 'end']
    num_extra = df.shape[1] - 3

    if extra_col_names:
        if len(extra_col_names) != num_extra:
            raise ValueError(f"Expected {num_extra} extra column names, got {len(extra_col_names)}")
        df.columns = base_cols + extra_col_names
    else:
        df.columns = base_cols + [f'col{i}' for i in range(4, 4 + num_extra)]

    return df


def expand_bed_regions(df: pd.DataFrame,
                       method: Optional[str] = None,
                       expand_bp: Optional[int] = None) -> pd.DataFrame:
    """
    Expand BED regions around center or edges.

    Parameters:
    - df (pd.DataFrame): BED-format DataFrame with columns 'chr', 'start', and 'end'.
    - method (str, optional): 'center' to expand symmetrically around midpoint,
                              'edge' to expand outward from start/end,
                              or None to return unmodified input.
    - expand_bp (int, optional): Number of base pairs to expand (required if method is specified).

    Returns:
    - pd.DataFrame: Modified DataFrame with adjusted 'start' and 'end' positions.

    Raises:
    - ValueError: If input BED is missing required columns.
    - ValueError: If method is specified but expand_bp is not.
    - ValueError: If method is not one of 'center', 'edge', or None.
    """
    df = df.copy()

    if not {"start", "end"}.issubset(df.columns):
        raise ValueError("DataFrame must contain 'start' and 'end' columns")

    if method is None:
        return df

    if expand_bp is None:
        raise ValueError("expand_bp must be provided if method is specified")

    if method == "center":
        midpoints = ((df["start"] + df["end"]) // 2).astype(int)
        df["start"] = (midpoints - expand_bp).clip(lower=0)
        df["end"] = midpoints + expand_bp
    elif method == "edge":
        df["start"] = (df["start"] - expand_bp).clip(lower=0)
        df["end"] = df["end"] + expand_bp
    else:
        raise ValueError("Method must be 'center', 'edge', or None")

    return df


def save_bed(df: pd.DataFrame, output_path: str) -> None:
    """
    Save a pandas DataFrame to a BED file.

    Parameters:
    - df (pd.DataFrame): DataFrame with at least 'chr', 'start', and 'end' columns.
    - output_path (str): Path to the output BED file.

    Raises:
    - ValueError: If required columns are missing.
    """
    required_cols = {"chr", "start", "end"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"DataFrame must contain columns: {required_cols}")

    df.to_csv(output_path, sep="\t", header=False, index=False)
