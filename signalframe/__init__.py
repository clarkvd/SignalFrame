import scipy
from packaging import version
import sys

required_scipy = "1.14.1"
if sys.version_info >= (3, 13) and version.parse(scipy.__version__) < version.parse(required_scipy):
    raise ImportError(
        f"signalframe requires scipy>={required_scipy} for Python 3.13+. "
        f"You have scipy {scipy.__version__}. Please upgrade using:\n\n"
        f"    pip install 'scipy=={required_scipy}'\n"
    )

from .io import load_bed, expand_bed_regions, save_bed
from .compute import compute_signal, compute_signal_multi
from .normalize import normalize_signal, compare_tracks, sort_signal_df
from .stats import compare_signal_groups, run_one_way_anova, run_two_way_anova
from .plot import plot_signal_region, plot_signals_from_bed

__all__ = [
    "load_bed", "expand_bed_regions", "save_bed",
    "compute_signal", "compute_signal_multi",
    "normalize_signal", "compare_tracks", "sort_signal_df",
    "compare_signal_groups", "run_one_way_anova", "run_two_way_anova",
    "plot_signal_region", "plot_signals_from_bed"
]
