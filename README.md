# SignalFrame

### `load_bed(bed_path: str, extra_col_names: Optional[List[str]] = None)`

**Description**  
Loads a BED file into a pandas DataFrame with standardized column names. Ensures the file contains at least three required columns: chromosome (`chr`), start (`start`), and end (`end`). If additional columns are present, they are automatically named `col4`, `col5`, etc., unless custom names are provided via the `extra_col_names` parameter.

**Parameters**
- `bed_path` (`str`): Path to the BED file.
- `extra_col_names` (`Optional[List[str]]`): Optional list of names for any columns beyond the first three.  
  Must match the number of extra columns if provided.

**Returns**
- `pd.DataFrame`: A DataFrame with at least the columns `chr`, `start`, and `end`.  
  Any additional columns are named automatically (`col4`, `col5`, ...) or using the provided names.

**Raises**
- `ValueError`: If the file has fewer than 3 columns.
- `ValueError`: If the number of `extra_col_names` does not match the number of extra columns.

**Examples**
```python
from signalframe import load_bed

# Load BED file with default column naming
bed_df = load_bed("regions.bed")
print(bed_df.head())

# Load BED file with custom column names
bed_df = load_bed("regions.bed", extra_col_names=["name", "score", "strand"])
print(bed_df.head())
```

### `expand_bed_regions(df: pd.DataFrame, method: Optional[str] = None, expand_bp: Optional[int] = None)`


**Description**  
Expands genomic intervals in a BED-format DataFrame either symmetrically around the center or outward from the edges. Useful for extending peak regions or extracting flanking sequence windows.

**Parameters**
- `df` (`pd.DataFrame`): A DataFrame containing at least the columns `'chr'`, `'start'`, and `'end'`.
- `method` (`str`, optional): Expansion strategy:
  - `'center'`: Expand symmetrically around the midpoint of each region.
  - `'edge'`: Expand outward from both the start and end coordinates.
  - `None`: Return the original DataFrame without modification.
- `expand_bp` (`int`, optional): Number of base pairs to expand on each side (required if `method` is specified).

**Returns**
- `pd.DataFrame`: A new DataFrame with adjusted `'start'` and `'end'` coordinates according to the specified expansion strategy.

**Raises**
- `ValueError`: If the input DataFrame is missing `'start'` or `'end'` columns.
- `ValueError`: If `method` is specified but `expand_bp` is not provided.
- `ValueError`: If `method` is not one of `'center'`, `'edge'`, or `None`.

**Examples**
```python
from signalframe import expand_bed_regions

# Example 1: Expand regions by 100 bp around the center
expanded_df = expand_bed_regions(bed_df, method="center", expand_bp=100)

# Example 2: Expand regions by 50 bp from the edges
expanded_df = expand_bed_regions(bed_df, method="edge", expand_bp=50)

# Example 3: Return unmodified regions
unchanged_df = expand_bed_regions(bed_df)
