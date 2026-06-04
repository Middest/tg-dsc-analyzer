"""Kinetics analysis: Coats-Redfern (single heating rate) and
iso-conversional methods (KAS, FWO, Friedman — multi heating rates).
"""

import numpy as np
import pandas as pd
from scipy.stats import linregress

R = 8.314  # J/(mol*K)


def coats_redfern(
    df: pd.DataFrame,
    T_range: tuple[float, float],
    alpha_range: tuple[float, float] = (0.05, 0.95),
    heating_rate_K_per_s: float | None = None,
    reaction_order: int = 1,
) -> dict:
    """Coats-Redfern first-order kinetic analysis for a single heating rate.

    For first-order: ln(-ln(1-alpha)/T^2) vs 1/T
    Ea = -slope * R

    Parameters
    ----------
    df : DataFrame with Temp_C, TG_percent columns.
    T_range : (T_low, T_high) range to analyze.
    alpha_range : fraction of conversion range to fit (avoids boundary instability).
    heating_rate_K_per_s : heating rate in K/s. If None, estimated from data.
    reaction_order : reaction order n (1 = first-order).

    Returns
    -------
    dict with Ea_kJ_mol, R2, slope, intercept, n_points used.
    """
    mask = (df["Temp_C"] >= T_range[0]) & (df["Temp_C"] <= T_range[1])
    sub = df.loc[mask].copy()
    if len(sub) < 10:
        return {"Ea_kJ_mol": np.nan, "R2": np.nan, "n_points": 0}

    T_C = sub["Temp_C"].values
    TG = sub["TG_percent"].values
    m0 = TG[0]
    m_final = TG[-1]

    if m0 - m_final <= 0:
        return {"Ea_kJ_mol": np.nan, "R2": np.nan, "n_points": len(T_C)}

    alpha = (m0 - TG) / (m0 - m_final)
    alpha = np.clip(alpha, 0.001, 0.999)

    a_mask = (alpha >= alpha_range[0]) & (alpha <= alpha_range[1])
    if a_mask.sum() < 5:
        # Fall back to wider range
        a_mask = (alpha >= 0.02) & (alpha <= 0.98)
        if a_mask.sum() < 5:
            return {"Ea_kJ_mol": np.nan, "R2": np.nan, "n_points": int(a_mask.sum())}

    alpha_fit = alpha[a_mask]
    T_K_fit = T_C[a_mask] + 273.15

    # Coats-Redfern for 1st order:
    # ln(-ln(1-alpha)/T^2) = ln(AR/bE) - E/(RT)
    # Y = ln(-ln(1-alpha)/T^2), X = 1/T
    # Slope = -E/R -> E = -slope * R
    # For nth order (n != 1):
    # ln((1-(1-alpha)^(1-n))/((1-n)*T^2)) = ln(AR/bE) - E/(RT)
    if reaction_order == 1:
        Y = np.log(-np.log(1 - alpha_fit) / (T_K_fit ** 2))
    else:
        n = reaction_order
        Y = np.log((1 - (1 - alpha_fit) ** (1 - n)) / ((1 - n) * T_K_fit ** 2))

    X = 1.0 / T_K_fit

    # Remove inf/nan
    valid = np.isfinite(Y) & np.isfinite(X)
    if valid.sum() < 5:
        return {"Ea_kJ_mol": np.nan, "R2": np.nan, "n_points": int(valid.sum())}

    slope, intercept, r_value, p_value, std_err = linregress(X[valid], Y[valid])

    Ea = -slope * R / 1000  # J/mol -> kJ/mol

    return {
        "Ea_kJ_mol": float(Ea),
        "R2": float(r_value ** 2),
        "slope": float(slope),
        "intercept": float(intercept),
        "n_points": int(valid.sum()),
        "T_range_K": (T_range[0] + 273.15, T_range[1] + 273.15),
    }


def pool_kinetics(
    df: pd.DataFrame,
    pools_scheme: str | list[dict] = "filimonenko2025",
    alpha_range: tuple[float, float] = (0.05, 0.95),
) -> pd.DataFrame:
    """Run Coats-Redfern kinetics for each defined thermal pool."""
    from .pools import get_pool_scheme
    pools = get_pool_scheme(pools_scheme)

    results = []
    for p in pools:
        cr = coats_redfern(df, (p["T_low"], p["T_high"]), alpha_range)
        cr["pool"] = p["pool"]
        cr["T_low_C"] = p["T_low"]
        cr["T_high_C"] = p["T_high"]
        results.append(cr)

    return pd.DataFrame(results)


def weighted_ea(
    pool_kinetics_df: pd.DataFrame,
    pool_mass_losses: dict[str, float],
) -> float:
    """Calculate weighted average Ea by pool mass-loss proportion.

    Parameters
    ----------
    pool_kinetics_df : DataFrame from pool_kinetics with 'pool' and 'Ea_kJ_mol' columns.
    pool_mass_losses : dict {pool_name: mass_loss_in_g_kg}.

    Returns
    -------
    float, weighted Ea in kJ/mol.
    """
    total_mass = sum(pool_mass_losses.values())
    if total_mass <= 0:
        return np.nan

    weighted = 0.0
    for _, row in pool_kinetics_df.iterrows():
        pool = row["pool"]
        ea = row["Ea_kJ_mol"]
        mass = pool_mass_losses.get(pool, 0)
        if not np.isnan(ea):
            weighted += ea * mass

    return weighted / total_mass


# --- iso-conversional methods (multi heating rate) ---

def _prepare_iso_data(
    samples: dict[str, pd.DataFrame],
    heating_rates: dict[str, float],
    T_range: tuple[float, float],
    n_T: int = 200,
) -> dict:
    """Prepare conversion-vs-temperature data across heating rates."""
    import warnings

    T_common = None
    alpha_dict = {}

    for name, df in samples.items():
        if name not in heating_rates:
            continue
        mask = (df["Temp_C"] >= T_range[0]) & (df["Temp_C"] <= T_range[1])
        sub = df.loc[mask]
        TG = sub["TG_percent"].values
        m0 = TG[0]
        m_final = TG[-1]
        alpha = (m0 - TG) / (m0 - m_final)
        alpha = np.clip(alpha, 0.001, 0.999)

        T = sub["Temp_C"].values
        if T_common is None:
            T_common = np.linspace(T.min(), T.max(), n_T)

        from scipy.interpolate import interp1d
        f = interp1d(
            alpha, T, kind="linear",
            bounds_error=False, fill_value="extrapolate",
        )
        T_at_alpha = f(np.linspace(0.05, 0.95, n_T))
        alpha_dict[name] = {
            "T": T_at_alpha,
            "beta": heating_rates[name],
            "alpha_vals": np.linspace(0.05, 0.95, n_T),
        }

    return alpha_dict


def kas(
    samples: dict[str, pd.DataFrame],
    heating_rates: dict[str, float],
    T_range: tuple[float, float],
) -> pd.DataFrame:
    """Kissinger-Akahira-Sunose (KAS) iso-conversional method.

    ln(beta/T^2) = const - Ea/(R*T)

    Requires multiple heating rates for the same material.
    """
    data = _prepare_iso_data(samples, heating_rates, T_range)

    n_points = len(next(iter(data.values()))["alpha_vals"])
    results = []

    for i in range(n_points):
        alpha = data[next(iter(data))]["alpha_vals"][i]
        X, Y = [], []
        for name, d in data.items():
            T_K = d["T"][i] + 273.15
            beta = d["beta"]  # K/min
            Y.append(np.log(beta / (T_K ** 2)))
            X.append(1.0 / T_K)

        if len(X) >= 3:
            slope, intercept, rv, _, _ = linregress(X, Y)
            Ea = -slope * R / 1000
            results.append({
                "alpha": float(alpha),
                "Ea_kJ_mol": float(Ea),
                "R2": float(rv ** 2),
                "T_mean_K": float(np.mean([d["T"][i] + 273.15 for d in data.values()])),
            })

    return pd.DataFrame(results)


def fwo(
    samples: dict[str, pd.DataFrame],
    heating_rates: dict[str, float],
    T_range: tuple[float, float],
) -> pd.DataFrame:
    """Flynn-Wall-Ozawa (FWO) iso-conversional method.

    ln(beta) = const - 1.052*Ea/(R*T)
    """
    data = _prepare_iso_data(samples, heating_rates, T_range)

    n_points = len(next(iter(data.values()))["alpha_vals"])
    results = []

    for i in range(n_points):
        alpha = data[next(iter(data))]["alpha_vals"][i]
        X, Y = [], []
        for name, d in data.items():
            T_K = d["T"][i] + 273.15
            beta = d["beta"]
            Y.append(np.log(beta))
            X.append(1.0 / T_K)

        if len(X) >= 3:
            slope, intercept, rv, _, _ = linregress(X, Y)
            Ea = -slope * R / 1.052 / 1000
            results.append({
                "alpha": float(alpha),
                "Ea_kJ_mol": float(Ea),
                "R2": float(rv ** 2),
            })

    return pd.DataFrame(results)
