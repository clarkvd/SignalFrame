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
import signalframe as sf

# Load BED file with default column naming
bed_df = sf.load_bed("regions.bed")
print(bed_df.head())

# Load BED file with custom column names
bed_df = sf.load_bed("regions.bed", extra_col_names=["name", "score", "strand"])
print(bed_df.head())
```

### `save_bed(df: pd.DataFrame, output_path: str)`

**Description**  
Saves a pandas DataFrame to a BED-format file using tab-separated values, without writing headers or index columns.

**Parameters**
- `df` (`pd.DataFrame`): DataFrame to save. Should ideally include at least `chr`, `start`, and `end` columns.
- `output_path` (`str`): Path to the output BED file.

**Returns**
- `None`

**Examples**
```python
import signalframe as sf

# Save DataFrame as BED file
sf.save_bed(bed_df, "output.bed")
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
import signalframe as sf

# Example 1: Expand regions by 100 bp around the center
expanded_df = sf.expand_bed_regions(bed_df, method="center", expand_bp=100)

# Example 2: Expand regions by 50 bp from the edges
expanded_df = sf.expand_bed_regions(bed_df, method="edge", expand_bp=50)

# Example 3: Return unmodified regions
unchanged_df = sf.expand_bed_regions(bed_df)
```

### `compute_signal(bigwig_path: str, bed_df: pd.DataFrame, method: str)`

**Description**  
Computes signal values from a BigWig file over a set of genomic intervals provided in BED format.  
The function supports multiple statistical methods (e.g., AUC, mean, max) and includes internal optimization via dynamic merging of nearby regions to reduce I/O operations and accelerate runtime.

**Parameters**
- `bigwig_path` (`str`): Path to the BigWig (.bw) file containing the signal track.
- `bed_df` (`pd.DataFrame`): DataFrame containing BED-format genomic intervals.  
  Must include `'chr'`, `'start'`, and `'end'` columns.
- `method` (`str`): Type of signal computation to perform. Must be one of:
  - `'auc'`: Sum of signal values across the region.
  - `'mean'`: Mean signal value.
  - `'max'`: Maximum signal value.
  - `'min'`: Minimum signal value.
  - `'median'`: Median signal value.
  - `'std'`: Standard deviation of signal values.
  - `'coverage'`: Count of non-zero signal positions.
  - `'nonzero_mean'`: Mean of non-zero signal values only.

**Returns**
- `pd.DataFrame`: The original DataFrame with an additional column named after the chosen `method` (e.g., `"auc"` or `"mean"`), containing the computed signal values for each region.

**Raises**
- `ValueError`: If the input DataFrame does not contain the required columns (`'chr'`, `'start'`, `'end'`).
- `ValueError`: If the specified `method` is not supported.

**Examples**
```python
import signalframe as sf

# Load BED DataFrame
bed_df = sf.load_bed("regions.bed")

# Compute area under the curve (AUC) for each region
bed_df_with_auc = sf.compute_signal("H3K27ac.bw", bed_df, method="auc")

# Compute mean signal
bed_df_with_mean = sf.compute_signal("ATAC_signal.bw", bed_df, method="mean")
```

### `compute_signal_multi(bigwig_paths: List[str], bed_df: pd.DataFrame, method: str)`

**Description**  
Computes signal values across multiple BigWig files for a common set of genomic intervals.  
Each BigWig track is processed independently, and the specified summary statistic is computed for every region in the BED file.  
The function merges nearby intervals internally to reduce I/O and optimize performance.

**Parameters**
- `bigwig_paths` (`List[str]`): List of paths to BigWig (.bw) signal files.
- `bed_df` (`pd.DataFrame`): DataFrame containing genomic intervals. Must include `'chr'`, `'start'`, and `'end'` columns.
- `method` (`str`): Signal summary statistic to compute. One of:
  - `'auc'`: Sum of signal values across the region.
  - `'mean'`: Mean signal value.
  - `'max'`: Maximum signal value.
  - `'min'`: Minimum signal value.
  - `'median'`: Median signal value.
  - `'std'`: Standard deviation of signal values.
  - `'coverage'`: Number of non-zero positions.
  - `'nonzero_mean'`: Mean of non-zero signal values only.

**Returns**
- `pd.DataFrame`: The original BED DataFrame with one new column per BigWig file, named as  
  `<basename>_<method>` (e.g., `H3K27ac_auc`, `ATAC_mean`, etc.).

**Raises**
- `ValueError`: If `method` is not valid or if required columns are missing.
- `ValueError`: If any of the provided BigWig files is invalid.

**Examples**
```python
import signalframe as sf

# Load a BED file
bed_df = sf.load_bed("peaks.bed")

# Define BigWig files
bw_files = ["H3K27ac.bw", "ATAC.bw", "H3K4me1.bw"]

# Compute AUC for each region in each file
result_df = sf.compute_signal_multi(bw_files, bed_df, method="auc")

# Compute mean signal
result_df = sf.compute_signal_multi(bw_files, bed_df, method="mean")
```

### `normalize_signal(df: pd.DataFrame, columns: Union[str, List[str]], method: str, pseudocount: float = 0.1, reference_matrix: Optional[Union[str, np.ndarray]] = None)`

**Description**  
Normalizes one or more signal columns in a DataFrame using common bioinformatics transformations.  
This is useful for standardizing signal values (e.g., AUCs) across different experiments, tracks, or conditions.

**Parameters**
- `df` (`pd.DataFrame`): DataFrame containing signal values to normalize.
- `columns` (`str` or `List[str]`): Name(s) of the column(s) to normalize.
- `method` (`str`): Normalization method. One of:
  - `'length'`: Normalize by region length (AUC per base pair).
  - `'zscore'`: Standard Z-score normalization.
  - `'log2'`: Log2 transform using `log2(x + pseudocount)`.
  - `'minmax'`: Min-max scale to [0, 1].
  - `'quantile'`: Quantile normalize to a reference distribution.
    - If `reference_matrix == "median"`, use the median across columns.
    - If `reference_matrix == "mean"`, use the mean across columns.
    - If `reference_matrix` is an `np.ndarray`, shape must match `(n_rows, n_columns)`.

- `pseudocount` (`float`, default = 0.1): Small constant added to avoid log of zero in `log2` normalization.
- `reference_matrix` (`str` or `np.ndarray`, optional): Used only with `'quantile'` normalization.

**Returns**
- `pd.DataFrame`: A new DataFrame containing the original columns and new columns named `<column>_norm` for each normalized signal column.

**Raises**
- `ValueError`: If column names are invalid, required columns are missing, or method/shape constraints are violated.

**Examples**
```python
import signalframe as sf

# Normalize by region length (AUC per bp)
df = sf.normalize_signal(df, columns="H3K27ac_auc", method="length")

# Z-score normalization
df = sf.normalize_signal(df, columns=["H3K27ac_auc", "ATAC_auc"], method="zscore")

# Log2 normalization with pseudocount
df = sf.normalize_signal(df, columns="H3K27ac_auc", method="log2", pseudocount=0.01)

# Min-max scaling
df = sf.normalize_signal(df, columns="ATAC_auc", method="minmax")

# Quantile normalization using median profile
df = sf.normalize_signal(df, columns=["H3K27ac_auc", "ATAC_auc"], method="quantile", reference_matrix="median")
```

### `compare_tracks(df: pd.DataFrame, reference: str, comparisons: Union[str, List[str]], mode: Union[str, List[str]] = ["difference", "fold_change", "log2FC", "percent_change"], pseudocount: float = 0.1)`

**Description**  
Compares one reference signal track to one or more other signal tracks across genomic regions.  
Generates new columns based on selected comparison metrics (e.g., fold change, log2 fold change, percent change).  
This is useful for contrasting peak intensity or accessibility between experimental conditions.

**Parameters**
- `df` (`pd.DataFrame`): DataFrame containing the signal values to compare.
- `reference` (`str`): Name of the reference signal column.
- `comparisons` (`str` or `List[str]`): Name(s) of target columns to compare against the reference.
- `mode` (`str` or `List[str]`): One or more comparison modes to compute:
  - `'difference'`: Absolute difference (`reference - target`)
  - `'fold_change'`: Fold change = `(reference + pseudocount) / (target + pseudocount)`
  - `'log2FC'`: Log2 fold change = `log2((reference + pseudocount) / (target + pseudocount))`
  - `'percent_change'`: Percent change = `(reference - target) / (target + pseudocount)`
- `pseudocount` (`float`, default = 0.1): Value added to avoid division by zero and undefined log operations.

**Returns**
- `pd.DataFrame`: A new DataFrame with columns added for each requested comparison.

**Raises**
- `ValueError`: If the reference or comparison columns are not found, or an invalid comparison mode is specified.

**Examples**
```python
import signalframe as sf

# Compare ATAC signal to H3K27ac using fold change and log2FC
df = sf.compare_tracks(df, reference="ATAC_auc", comparisons="H3K27ac_auc", mode=["fold_change", "log2FC"])

# Compare against multiple tracks using all metrics
df = sf.compare_tracks(df, reference="ATAC_auc", comparisons=["H3K27ac_auc", "H3K4me1_auc"])
```

### `sort_signal_df(df: pd.DataFrame, sort_by: str = "genomic_position", ascending: bool = True)`

**Description**  
Sorts a DataFrame containing genomic signal values either by genomic position (`chr`, `start`) or by any specified column.  
This is useful for organizing results before visualization, exporting, or downstream analysis.

**Parameters**
- `df` (`pd.DataFrame`): DataFrame containing signal values and genomic intervals.
- `sort_by` (`str`, default = `"genomic_position"`): Sorting mode:
  - `"genomic_position"`: Sorts by `'chr'` and `'start'` columns.
  - Any other string is interpreted as a column name to sort by.
- `ascending` (`bool`, default = `True`): Sort order.
  - `True`: Sort in ascending order.
  - `False`: Sort in descending order.

**Returns**
- `pd.DataFrame`: A sorted copy of the input DataFrame with index reset.

**Raises**
- `ValueError`: If `sort_by` is not `"genomic_position"` and is not a valid column in the DataFrame.
- `ValueError`: If `"genomic_position"` is chosen but `'chr'` and `'start'` columns are missing.

**Examples**
```python
import signalframe as sf

# Sort by genomic coordinates
sorted_df = sf.sort_signal_df(df, sort_by="genomic_position")

# Sort by signal intensity (e.g., AUC column)
sorted_df = sf.sort_signal_df(df, sort_by="H3K27ac_auc", ascending=False)
```


### `compare_signal_groups(df: pd.DataFrame, group_col: str, value_col: str, test: Literal["t-test", "mannwhitney"] = "t-test")`

**Description**  
Performs a statistical comparison of signal values between two groups using either an unpaired t-test or the non-parametric Mann–Whitney U test.  
Useful for evaluating differential signals (e.g., peak intensities) between experimental conditions.

**Parameters**
- `df` (`pd.DataFrame`): DataFrame containing the data to compare.
- `group_col` (`str`): Column name indicating the group for each row. Must contain exactly two unique values.
- `value_col` (`str`): Column containing the numeric signal values to be compared.
- `test` (`str`, default = `"t-test"`): Statistical test to perform. Options:
  - `"t-test"`: Welch’s t-test (assumes unequal variances).
  - `"mannwhitney"`: Mann–Whitney U test (non-parametric).

**Returns**
- `float`: The p-value resulting from the statistical test.

**Raises**
- `ValueError`: If `group_col` does not contain exactly two unique groups.
- `ValueError`: If `test` is not one of `"t-test"` or `"mannwhitney"`.

**Examples**
```python
import signalframe as sf

# Compare mean signal between FoxP3 and input groups using t-test
pval = sf.compare_signal_groups(df, group_col="group", value_col="H3K27ac_auc", test="t-test")

# Use Mann–Whitney U test instead
pval = sf.compare_signal_groups(df, group_col="group", value_col="ATAC_auc", test="mannwhitney")
```


### `run_one_way_anova(df: pd.DataFrame, factor: str, response: str, return_model: bool = False)`

**Description**  
Performs a one-way ANOVA (Analysis of Variance) to test whether the means of a numeric response differ significantly across groups defined by a single categorical factor.

Uses `statsmodels` under the hood to fit an OLS model and generate the ANOVA summary table.

**Parameters**
- `df` (`pd.DataFrame`): Input DataFrame containing the data.
- `factor` (`str`): Name of the column representing the categorical grouping variable.
- `response` (`str`): Name of the column representing the numeric response variable.
- `return_model` (`bool`, default = `False`): If `True`, also return the fitted statsmodels OLS model.

**Returns**
- `pd.DataFrame`: ANOVA summary table (F-statistic, p-value, degrees of freedom, etc.).
- `RegressionResultsWrapper` (optional): The fitted OLS model if `return_model=True`.

**Raises**
- `ValueError`: If either the `factor` or `response` column is missing from the DataFrame.

**Examples**
```python
import signalframe as sf

# Run one-way ANOVA comparing H3K27ac signal across tissue types
anova_table = sf.run_one_way_anova(df, factor="tissue", response="H3K27ac_auc")

# Also retrieve the fitted model for further inspection
anova_table, model = sf.run_one_way_anova(df, factor="condition", response="ATAC_signal", return_model=True)
```


### `run_two_way_anova(df: pd.DataFrame, factor1: str, factor2: str, response: str, include_interaction: bool = True, return_model: bool = False)`

**Description**  
Performs a two-way ANOVA (Analysis of Variance) to test the effects of two categorical factors on a numeric response variable.  
Optionally includes an interaction term to evaluate whether the effect of one factor depends on the level of the other.

Uses `statsmodels` OLS modeling and ANOVA Type II sum of squares.

**Parameters**
- `df` (`pd.DataFrame`): Input DataFrame with the data.
- `factor1` (`str`): First categorical factor column name.
- `factor2` (`str`): Second categorical factor column name.
- `response` (`str`): Name of the numeric response column.
- `include_interaction` (`bool`, default = `True`): Whether to include the interaction term `C(factor1):C(factor2)` in the model.
- `return_model` (`bool`, default = `False`): If `True`, return both the ANOVA table and the fitted model object.

**Returns**
- `pd.DataFrame`: ANOVA summary table showing effects of each factor and (optionally) their interaction.
- `RegressionResultsWrapper` (optional): The fitted OLS model if `return_model=True`.

**Raises**
- `ValueError`: If any of the specified columns are missing from the DataFrame.

**Examples**
```python
import signalframe as sf

# Run two-way ANOVA including interaction
anova_table = sf.run_two_way_anova(df, factor1="tissue", factor2="condition", response="H3K27ac_auc")

# Run without interaction and return the fitted model as well
anova_table, model = sf.run_two_way_anova(
    df,
    factor1="cell_type",
    factor2="treatment",
    response="signal",
    include_interaction=False,
    return_model=True
)
```

### `plot_signal_region(bigwig_paths: Union[str, List[str]], chrom: str, start: int, end: int, labels: Optional[Union[str, List[str]]] = None, title: Optional[str] = None, figsize: tuple = (12, 4))`

**Description**  
Plots the signal values from one or more BigWig files over a specified genomic region.  
Useful for visual inspection of peak shapes, comparisons between tracks, and identifying signal differences across datasets.

**Parameters**
- `bigwig_paths` (`str` or `List[str]`): Path(s) to BigWig file(s).
- `chrom` (`str`): Chromosome name (e.g., `"chr1"`).
- `start` (`int`): Start coordinate (0-based).
- `end` (`int`): End coordinate (non-inclusive).
- `labels` (`str` or `List[str]`, optional): Optional labels for each track.  
  Defaults to the file name(s) if not provided.
- `title` (`str`, optional): Title for the plot.
- `figsize` (`tuple`, default = `(12, 4)`): Size of the output figure in inches.

**Raises**
- `ValueError`: If the number of labels does not match the number of BigWig files.
- Prints warnings if chromosomes are not found or files cannot be read.

**Returns**
- Displays a matplotlib plot of signal values across the specified region.

**Examples**
```python
import signalframe as sf

# Plot a single BigWig track
sf.plot_signal_region("ATAC.bw", chrom="chr15", start=100000, end=101000)

# Plot multiple tracks with custom labels and a title
sf.plot_signal_region(
    bigwig_paths=["H3K27ac.bw", "ATAC.bw"],
    chrom="chr7",
    start=50000,
    end=52000,
    labels=["H3K27ac", "ATAC-seq"],
    title="Signal at chr7:50,000-52,000"
)
```


### `plot_signals_from_bed(bed_df: pd.DataFrame, bigwig_paths: Union[str, List[str]], max_plots: int = 10, shared_y: bool = True, labels: Optional[Union[str, List[str]]] = None, figsize: tuple = (12, 3))`

**Description**  
Generates a series of signal plots from one or more BigWig files across genomic regions defined in a BED-format DataFrame.  
This function is ideal for visually inspecting signal profiles over peaks, motifs, or other custom regions.

**Parameters**
- `bed_df` (`pd.DataFrame`): DataFrame with BED-style columns `'chr'`, `'start'`, and `'end'`.
- `bigwig_paths` (`str` or `List[str]`): Path(s) to one or more BigWig signal files.
- `max_plots` (`int`, default = `10`): Maximum number of regions (rows) to visualize.
- `shared_y` (`bool`, default = `True`): Whether to use the same y-axis scale across all subplots.
- `labels` (`str` or `List[str]`, optional): Optional custom labels for each BigWig track.  
  Defaults to the base filenames if not specified.
- `figsize` (`tuple`, default = `(12, 3)`): Size of each individual subplot (width, height in inches).

**Returns**
- Displays a series of stacked line plots, one for each region (up to `max_plots`).

**Raises**
- `ValueError`: If the number of labels does not match the number of BigWig files.

**Notes**
- If `shared_y=True`, the function first computes global y-axis limits across all selected regions.
- Regions with missing chromosomes in any BigWig file are skipped with a warning.

**Examples**
```python
import signalframe as sf

# Plot the first 10 peaks in the BED DataFrame
sf.plot_signals_from_bed(bed_df, bigwig_paths="H3K27ac.bw")

# Plot up to 5 regions with two BigWig tracks, shared y-axis turned off
sf.plot_signals_from_bed(
    bed_df,
    bigwig_paths=["ATAC.bw", "H3K27ac.bw"],
    max_plots=5,
    shared_y=False,
    labels=["ATAC-seq", "H3K27ac"]
)
```

