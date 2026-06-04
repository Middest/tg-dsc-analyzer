"""Preprocessing: smoothing, baseline correction, interpolation, deduplication."""

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d


def smooth_savgol(
    df: pd.DataFrame,
    columns: list[str] | None = None,
    window: int = 15,
    polyorder: int = 3,
) -> pd.DataFrame:
    """Apply Savitzky-Golay smoothing to specified columns.

    Parameters
    ----------
    df : DataFrame with Temp_C and signal columns.
    columns : list of column names to smooth. Defaults to TG, DTG, DSC.
    window : odd integer, smoothing window size.
    polyorder : polynomial order.
    """
    if columns is None:
        columns = ["TG_percent", "DTG_percent_per_C", "DSC_uW_per_mg"]
    if window % 2 == 0:
        window += 1
    if window < polyorder + 2:
        window = polyorder + 2
        if window % 2 == 0:
            window += 1

    out = df.copy()
    for col in columns:
        if col in out.columns:
            out[col] = savgol_filter(out[col].values, window, polyorder)
    return out


def baseline_correct_dsc(
    df: pd.DataFrame,
    T_range: tuple[float, float] | None = None,
    method: str = "linear",
) -> pd.DataFrame:
    """Perform DSC baseline correction.

    Parameters
    ----------
    df : DataFrame with Temp_C and DSC_uW_per_mg.
    T_range : (T_low, T_high) region to use for baseline anchor points.
    method : "linear" only for now.
    """
    out = df.copy()
    if T_range is None:
        mask = np.ones(len(df), dtype=bool)
    else:
        mask = (df["Temp_C"] >= T_range[0]) & (df["Temp_C"] <= T_range[1])

    T = df["Temp_C"].values[mask]
    dsc = df["DSC_uW_per_mg"].values[mask]

    if method == "linear":
        # Linear baseline between endpoints
        baseline = np.interp(
            df["Temp_C"].values, [T[0], T[-1]], [dsc[0], dsc[-1]]
        )
        out["DSC_uW_per_mg"] = df["DSC_uW_per_mg"].values - baseline

    return out


def interpolate_common_temp(
    samples: dict[str, pd.DataFrame],
    T_range: tuple[float, float] | None = None,
    n_points: int = 1024,
) -> dict[str, pd.DataFrame]:
    """Interpolate all samples onto a common temperature grid."""
    if T_range is None:
        T_min = max(s["Temp_C"].min() for s in samples.values())
        T_max = min(s["Temp_C"].max() for s in samples.values())
    else:
        T_min, T_max = T_range

    T_common = np.linspace(T_min, T_max, n_points)

    result = {}
    for name, df in samples.items():
        interp_df = pd.DataFrame({"Temp_C": T_common})
        for col in ["TG_percent", "DTG_percent_per_C", "DSC_uW_per_mg"]:
            if col in df.columns:
                f = interp1d(
                    df["Temp_C"].values, df[col].values,
                    kind="linear", bounds_error=False,
                    fill_value=(df[col].values[0], df[col].values[-1]),
                )
                interp_df[col] = f(T_common)
        interp_df["Time_min"] = np.linspace(
            df["Time_min"].min(), df["Time_min"].max(), n_points
        )
        for col in df.columns:
            if col not in interp_df.columns:
                interp_df[col] = df[col].iloc[0] if len(df) > 0 else None
        result[name] = interp_df

    return result


def remove_duplicate_time(df: pd.DataFrame) -> pd.DataFrame:
    """Remove rows with duplicate time values (keep first)."""
    if "Time_min" in df.columns:
        return df.drop_duplicates(subset="Time_min", keep="first").reset_index(drop=True)
    return df
