"""DSC analysis: heat flow curves, peak detection, and energy integration."""

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.signal import find_peaks, peak_widths


def heat_flow(df: pd.DataFrame) -> np.ndarray:
    """Return DSC heat flow signal (uW/mg)."""
    return df["DSC_uW_per_mg"].values


def detect_peaks(
    df: pd.DataFrame,
    T_range: tuple[float, float] | None = None,
    height: float | None = None,
    prominence: float | None = None,
    direction: str = "both",
) -> list[dict]:
    """Detect endothermic/exothermic peaks in DSC signal.

    Returns list of dicts with peak properties.
    """
    temp = df["Temp_C"].values
    dsc = df["DSC_uW_per_mg"].values

    if T_range is not None:
        mask = (temp >= T_range[0]) & (temp <= T_range[1])
        temp = temp[mask]
        dsc = dsc[mask]

    if height is None:
        height = 0.05 * (np.max(dsc) - np.min(dsc))
    if prominence is None:
        prominence = 0.02 * (np.max(dsc) - np.min(dsc))

    peaks = []
    if direction in ("both", "exo"):
        # Exothermic = positive (for NETZSCH convention, DSC > 0)
        pos_peaks, props = find_peaks(dsc, height=height, prominence=prominence)
        widths = peak_widths(dsc, pos_peaks, rel_height=0.5)
        for i, idx in enumerate(pos_peaks):
            peaks.append({
                "T_peak_C": float(temp[idx]),
                "DSC_peak_uW_mg": float(dsc[idx]),
                "type": "exothermic",
                "width_C": float(widths[0][i]) if i < len(widths[0]) else None,
                "T_start_C": float(temp[int(widths[2][i])]) if i < len(widths[2]) else None,
                "T_end_C": float(temp[int(widths[3][i])]) if i < len(widths[3]) else None,
            })

    if direction in ("both", "endo"):
        neg_peaks, props = find_peaks(-dsc, height=height, prominence=prominence)
        widths = peak_widths(-dsc, neg_peaks, rel_height=0.5)
        for i, idx in enumerate(neg_peaks):
            peaks.append({
                "T_peak_C": float(temp[idx]),
                "DSC_peak_uW_mg": float(dsc[idx]),
                "type": "endothermic",
                "width_C": float(widths[0][i]) if i < len(widths[0]) else None,
                "T_start_C": float(temp[int(widths[2][i])]) if i < len(widths[2]) else None,
                "T_end_C": float(temp[int(widths[3][i])]) if i < len(widths[3]) else None,
            })

    return sorted(peaks, key=lambda p: p["T_peak_C"])


def peak_integral(
    df: pd.DataFrame,
    T_start: float,
    T_end: float,
    signal_col: str = "DSC_uW_per_mg",
) -> dict:
    """Integrate DSC signal over a temperature range.

    DSC is in uW/mg = uJ/(s*mg). Integration over time (seconds) gives uJ/mg.
    """
    mask = (df["Temp_C"] >= T_start) & (df["Temp_C"] <= T_end)
    if mask.sum() < 2:
        return {"energy_uJ_mg": np.nan, "T_start": T_start, "T_end": T_end, "n_points": int(mask.sum())}

    time_s = df.loc[mask, "Time_min"].values * 60  # min -> s
    signal = df.loc[mask, signal_col].values

    energy = trapezoid(signal, time_s)
    return {
        "energy_uJ_mg": float(energy),
        "energy_J_g_soil": float(energy * 1e-3),  # uJ/mg -> J/g
        "T_start": T_start,
        "T_end": T_end,
        "n_points": int(mask.sum()),
    }


def total_energy(
    df: pd.DataFrame,
    T_range: tuple[float, float],
    signal_col: str = "DSC_uW_per_mg",
) -> dict:
    """Total energy integral over a temperature range."""
    return peak_integral(df, T_range[0], T_range[1], signal_col)


def energy_density(
    df: pd.DataFrame,
    T_range: tuple[float, float],
) -> dict:
    """Energy density ED = DSC integral / TG mass loss (kJ/g OM).

    This is a key parameter indicating the energy stored per unit organic matter lost.

    Physical plausibility checks:
      - ED should typically be 10–500 kJ/g OM for soil/biochar samples
      - Negative ED indicates DSC baseline drift or integration error
      - ED > 1000 kJ/g OM suggests near-zero mass loss (division instability)
    """
    from .tg import mass_loss as tg_mass_loss

    energy = total_energy(df, T_range)["energy_J_g_soil"]
    mass_loss_pct = tg_mass_loss(df, T_range[0], T_range[1])
    mass_loss_frac = mass_loss_pct / 100  # percent -> fraction

    warnings = []

    if mass_loss_frac <= 0:
        warnings.append("mass_loss_zero")
        return {
            "ED_kJ_g_OM": np.nan, "energy_J_g_soil": energy,
            "mass_loss_pct": mass_loss_pct,
            "ed_warnings": warnings,
        }

    # ED_kJ_g_OM = E(J/g_soil) / mass_loss_fraction(g_OM/g_soil) / 1000
    ed = energy / mass_loss_frac / 1000

    if ed < 0:
        warnings.append("ED_negative")
    if abs(ed) > 1000:
        warnings.append("ED_extreme")
    if mass_loss_pct < 0.5:
        warnings.append("mass_loss_very_low")

    return {
        "ED_kJ_g_OM": float(ed) if ed >= 0 else np.nan,
        "ED_kJ_g_OM_raw": float(ed),  # keep raw value for debugging
        "energy_J_g_soil": energy,
        "mass_loss_pct": mass_loss_pct,
        "T_range": T_range,
        "ed_warnings": ";".join(warnings) if warnings else "ok",
    }
