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
