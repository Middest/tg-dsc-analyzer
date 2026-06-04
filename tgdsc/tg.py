"""TG/DTG analysis: mass loss curves and characteristic temperatures."""

import numpy as np
import pandas as pd


def mass_loss(
    df: pd.DataFrame, T_start: float, T_end: float
) -> float:
    """Mass loss (percent) between two temperatures."""
    mask = (df["Temp_C"] >= T_start) & (df["Temp_C"] <= T_end)
    if mask.sum() < 2:
        return np.nan
    tg_vals = df.loc[mask, "TG_percent"].values
    return tg_vals[0] - tg_vals[-1]


def mass_loss_at_T(df: pd.DataFrame, T_target: float) -> float:
    """TG percent value at a specific temperature (interpolated)."""
    idx = np.searchsorted(df["Temp_C"].values, T_target)
    idx = np.clip(idx, 0, len(df) - 1)
    return df["TG_percent"].values[idx]


def cumulative_mass_loss_curve(
    df: pd.DataFrame, T_range: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray]:
    """Cumulative mass loss within a temperature range.

    Returns (Temp, cumulative_loss_pct) arrays.
    """
    mask = (df["Temp_C"] >= T_range[0]) & (df["Temp_C"] <= T_range[1])
    T = df.loc[mask, "Temp_C"].values
    tg = df.loc[mask, "TG_percent"].values
    cum_loss = tg[0] - tg
    return T, cum_loss


def T50(df: pd.DataFrame, T_range: tuple[float, float]) -> float | None:
    """Temperature at which 50% of cumulative mass loss within T_range is reached.

    This is the primary thermal stability metric used in SOM studies.
    """
    T, cum_loss = cumulative_mass_loss_curve(df, T_range)
    total = cum_loss[-1] - cum_loss[0]
    if total <= 0:
        return None
    half = total / 2
    idx = np.searchsorted(cum_loss, half)
    if idx >= len(T):
        return T[-1]
    # Linear interpolation
    if idx == 0:
        return T[0]
    frac = (half - cum_loss[idx - 1]) / (cum_loss[idx] - cum_loss[idx - 1])
    return float(T[idx - 1] + frac * (T[idx] - T[idx - 1]))


def char_temperatures(df: pd.DataFrame) -> dict:
    """Extract characteristic temperatures from DTG curve.

    Returns dict with T_onset, T_peak_dtg, T_end.
    """
    dtg = df["DTG_percent_per_C"].values
    T = df["Temp_C"].values

    peak_idx = np.argmax(np.abs(dtg))
    T_peak = T[peak_idx]

    # Onset: intersection of baseline tangent with max-slope tangent
    # Simplified: temperature where DTG first exceeds 5% of peak value
    threshold = 0.05 * np.abs(dtg[peak_idx])
    onset_idx = np.where(np.abs(dtg[:peak_idx]) >= threshold)[0]
    T_onset = T[onset_idx[0]] if len(onset_idx) > 0 else T[0]

    # End: DTG returns near zero after peak
    end_candidates = np.where(
        (np.abs(dtg[peak_idx:]) < threshold) & (T[peak_idx:] > T_peak)
    )[0]
    T_end = T[peak_idx + end_candidates[0]] if len(end_candidates) > 0 else T[-1]

    return {"T_onset_C": float(T_onset), "T_peak_dtg_C": float(T_peak), "T_end_C": float(T_end)}


def residual_mass_percent(df: pd.DataFrame, T: float) -> float:
    """Residual mass percent at a given temperature."""
    idx = np.searchsorted(df["Temp_C"].values, T)
    idx = np.clip(idx, 0, len(df) - 1)
    return float(df["TG_percent"].values[idx])


def fractional_conversion(
    df: pd.DataFrame, T_range: tuple[float, float]
) -> pd.DataFrame:
    """Compute fractional conversion alpha for kinetics analysis.

    alpha = (m0 - m_T) / (m0 - m_final) within T_range.
    """
    mask = (df["Temp_C"] >= T_range[0]) & (df["Temp_C"] <= T_range[1])
    sub = df.loc[mask].copy()
    m0 = sub["TG_percent"].values[0]
    m_final = sub["TG_percent"].values[-1]
    sub["alpha"] = (m0 - sub["TG_percent"].values) / (m0 - m_final)
    sub["alpha"] = sub["alpha"].clip(0, 1)
    return sub
