import os
import pyBigWig
import pandas as pd
import numpy as np
from typing import List, Union, Optional, Dict, Tuple
import matplotlib.pyplot as plt

def plot_signal_region(bigwig_paths: Union[str, List[str]],
                       chrom: str,
                       start: int,
                       end: int,
                       labels: Optional[Union[str, List[str]]] = None,
                       title: Optional[str] = None,
                       figsize: tuple = (12, 4)):
    """
    Plot signal from one or more BigWig files across a specified genomic region.

    Parameters:
    - bigwig_paths (str or List[str]): Path(s) to BigWig file(s).
    - chrom (str): Chromosome name (e.g., "chr15").
    - start (int): Start position (0-based).
    - end (int): End position (non-inclusive).
    - labels (str or List[str], optional): Plot labels for each track. Defaults to file basenames.
    - title (str, optional): Optional title for the plot.
    - figsize (tuple): Size of the matplotlib figure in inches.

    Raises:
    - ValueError: If label count does not match number of BigWig files.
    """
    if isinstance(bigwig_paths, str):
        bigwig_paths = [bigwig_paths]
    if labels is None:
        labels = [os.path.basename(p) for p in bigwig_paths]
    elif isinstance(labels, str):
        labels = [labels]

    if len(labels) != len(bigwig_paths):
        raise ValueError("Number of labels must match number of BigWig files.")

    plt.figure(figsize=figsize)

    for i, bw_path in enumerate(bigwig_paths):
        label = labels[i]
        try:
            with pyBigWig.open(bw_path) as bw:
                if chrom not in bw.chroms():
                    print(f"[Warning] Chromosome '{chrom}' not found in {label}. Skipping.")
                    continue
                signal = bw.values(chrom, start, end, numpy=True)
                signal = np.nan_to_num(signal, nan=0.0)
                x = np.arange(start, end)
                plt.plot(x, signal, label=label)
        except Exception as e:
            print(f"[Error] Failed to read {bw_path}: {e}")

    plt.xlabel(f"Genomic Position on {chrom}")
    plt.ylabel("Signal")
    if title:
        plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_signals_from_bed(bed_df: pd.DataFrame,
                          bigwig_paths: Union[str, List[str]],
                          max_plots: int = 10,
                          shared_y: bool = True,
                          labels: Optional[Union[str, List[str]]] = None,
                          figsize: tuple = (12, 3)):
    """
    Plot signal tracks from one or more BigWig files across regions in a BED dataframe.

    Parameters:
    - bed_df (pd.DataFrame): DataFrame with columns 'chr', 'start', and 'end'.
    - bigwig_paths (str or List[str]): Path(s) to BigWig file(s).
    - max_plots (int): Maximum number of regions (rows) from the BED file to plot. Default is 10.
    - shared_y (bool): If True, all subplots will share the same y-axis scale.
    - labels (str or List[str], optional): Labels for each BigWig file; defaults to file basenames.
    - figsize (tuple): Size of each subplot (width, height in inches).

    Raises:
    - ValueError: If the number of labels does not match the number of BigWig files.
    """
    # Normalize inputs
    if isinstance(bigwig_paths, str):
        bigwig_paths = [bigwig_paths]
    if labels is None:
        labels = [os.path.basename(p) for p in bigwig_paths]
    elif isinstance(labels, str):
        labels = [labels]

    if len(labels) != len(bigwig_paths):
        raise ValueError("Length of labels must match number of BigWig files.")

    n_plots = min(len(bed_df), max_plots)
    fig, axs = plt.subplots(n_plots, 1, figsize=(figsize[0], figsize[1] * n_plots), sharey=shared_y)

    if n_plots == 1:
        axs = [axs]

    bws = [pyBigWig.open(p) for p in bigwig_paths]

    # Optional global y-axis calculation
    global_ymin, global_ymax = np.inf, -np.inf
    if shared_y:
        for i in range(n_plots):
            chrom, start, end = bed_df.iloc[i][['chr', 'start', 'end']]
            for bw in bws:
                if chrom not in bw.chroms():
                    continue
                values = bw.values(chrom, start, end, numpy=True)
                values = np.nan_to_num(values, nan=0.0)
                if values.size > 0:
                    global_ymin = min(global_ymin, np.min(values))
                    global_ymax = max(global_ymax, np.max(values))
        if not np.isfinite(global_ymin) or not np.isfinite(global_ymax):
            global_ymin, global_ymax = 0, 1  # default safe fallback

    for i in range(n_plots):
        chrom, start, end = bed_df.iloc[i][['chr', 'start', 'end']]
        ax = axs[i]
        x = np.arange(start, end)
        for j, bw in enumerate(bws):
            label = labels[j]
            if chrom not in bw.chroms():
                print(f"[Warning] {chrom} not found in {label}")
                continue
            values = bw.values(chrom, start, end, numpy=True)
            values = np.nan_to_num(values, nan=0.0)
            ax.plot(x, values, label=label)

        ax.set_title(f"{chrom}:{start}-{end}")
        ax.set_xlabel("Genomic Position")
        ax.set_ylabel("Signal")
        if shared_y:
            ax.set_ylim(global_ymin, global_ymax)

    for bw in bws:
        bw.close()

    axs[0].legend(loc="upper right")
    plt.tight_layout()
    plt.show()