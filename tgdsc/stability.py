"""Thermal stability indices for soil organic matter.

Metrics:
  - T50: temperature at 50% cumulative mass loss
  - R50: residual mass after oxidation
  - TI: thermostability index = mass_loss(350-550) / mass_loss(200-550)
  - ED: energy density (kJ/g OM)
  - Comprehensive multi-index report
"""

import numpy as np
import pandas as pd


def T50(
    df: pd.DataFrame,
    T_range: tuple[float, float] = (190, 640),
) -> dict:
    """Temperature at 50% cumulative mass loss within T_range."""
    from .tg import T50 as _t50_calc
    val = _t50_calc(df, T_range)
    return {"T50_C": val, "T_range": T_range}


def R50(df: pd.DataFrame, T_ref: float = 800) -> dict:
    """Residual mass percent at reference temperature."""
    from .tg import residual_mass_percent
    idx = np.searchsorted(df["Temp_C"].values, T_ref)
    idx = np.clip(idx, 0, len(df) - 1)
    # R50 = (residual at T) / (initial mass) * 100
    # Clip to first value > 99% to avoid moisture artifact
    tg_start = df["TG_percent"].values[0]
    tg_ref = df["TG_percent"].values[idx]
    return {
        "R50_percent": float(tg_ref),
        "T_ref_C": T_ref,
        "TG_initial_percent": float(tg_start),
    }


def thermostability_index(
    df: pd.DataFrame,
    T1: float = 200,
    T2: float = 350,
    T3: float = 550,
) -> dict:
    """Thermostability Index (TI).

    TI = mass_loss(350-550°C) / mass_loss(200-550°C)

    Higher TI indicates more recalcitrant (thermally stable) organic matter.
    Standard reference: Gregorich et al., Plante et al.
    """
    from .tg import mass_loss

    ml_200_350 = mass_loss(df, T1, T2)
    ml_350_550 = mass_loss(df, T2, T3)
    ml_200_550 = ml_200_350 + ml_350_550

    ti = ml_350_550 / ml_200_550 if ml_200_550 > 0 else np.nan

    return {
        "TI": float(ti) if not np.isnan(ti) else np.nan,
        "mass_loss_200_350": ml_200_350,
        "mass_loss_350_550": ml_350_550,
        "mass_loss_200_550": ml_200_550,
        "TG_at_200C": float(np.interp(T1, df["Temp_C"].values, df["TG_percent"].values)),
        "TG_at_350C": float(np.interp(T2, df["Temp_C"].values, df["TG_percent"].values)),
        "TG_at_550C": float(np.interp(T3, df["Temp_C"].values, df["TG_percent"].values)),
    }


def energy_density(
    df: pd.DataFrame,
    T_range: tuple[float, float] = (190, 640),
) -> dict:
    """Energy Density (ED) - DSC energy integral per unit mass loss."""
    from .dsc import energy_density as _ed_func
    return _ed_func(df, T_range)


def comprehensive_report(
    df: pd.DataFrame,
    pools_scheme: str = "filimonenko2025",
    T_range: tuple[float, float] = (190, 640),
    include_kinetics: bool = True,
) -> dict:
    """Complete thermal stability assessment for one sample.

    Returns a dict with all major indices and pool-level metrics.
    """
    from .pools import calculate_pools, get_pool_scheme
    from .tg import mass_loss

    report = {}
    report["T_range"] = T_range

    # T50
    report.update(T50(df, T_range))

    # Mass loss over range
    report["SOM_loss_percent"] = mass_loss(df, T_range[0], T_range[1])
    report["SOM_loss_g_kg_soil"] = report["SOM_loss_percent"] * 10  # % of soil -> g/kg

    # ED
    report.update(energy_density(df, T_range))

    # Thermostability Index
    report.update(thermostability_index(df))

    # Pool analysis
    pools_df = calculate_pools(df, pools_scheme, T_range)
    for _, row in pools_df.iterrows():
        pfx = row["pool"]
        report[f"{pfx}_mass_loss_percent"] = row["mass_loss_percent"]
        report[f"{pfx}_SOM_pool_g_kg"] = row["SOM_pool_g_kg_soil"]
        report[f"{pfx}_proportion_percent"] = row["pool_proportion_percent"]
        report[f"{pfx}_ED_kJ_g_OM"] = row["ED_kJ_g_OM"]

    # Kinetics
    if include_kinetics:
        from .kinetics import pool_kinetics, weighted_ea
        kin_df = pool_kinetics(df, pools_scheme)
        pool_masses = {}
        for _, row in pools_df.iterrows():
            pool_masses[row["pool"]] = row["SOM_pool_g_kg_soil"]
        report["Ea_weighted_kJ_mol"] = weighted_ea(kin_df, pool_masses)

        for _, row in kin_df.iterrows():
            report[f"Ea_{row['pool']}_kJ_mol"] = row["Ea_kJ_mol"]
            report[f"Ea_{row['pool']}_R2"] = row["R2"]

    return report


def batch_report(
    samples: dict[str, pd.DataFrame],
    treatment_map: dict[str, dict] | None = None,
    pools_scheme: str = "filimonenko2025",
    T_range: tuple[float, float] = (190, 640),
) -> pd.DataFrame:
    """Run comprehensive_report for all samples, return DataFrame."""
    rows = []
    for name, df in samples.items():
        try:
            r = comprehensive_report(df, pools_scheme, T_range)
            r["sample"] = name
            if treatment_map and name in treatment_map:
                for k, v in treatment_map[name].items():
                    r[k] = v
            rows.append(r)
        except Exception as e:
            print(f"Error processing {name}: {e}")

    return pd.DataFrame(rows)


def compare_treatments(
    report_df: pd.DataFrame,
    control_label: str = "CK",
    treatment_col: str = "Treatment",
) -> pd.DataFrame:
    """Calculate relative response ratios vs control for all numeric columns.

    If treatment_col not found, falls back to 'sample' column.
    If control not found, returns original DataFrame with a warning.
    """
    # Fallback: use 'sample' as treatment column
    if treatment_col not in report_df.columns:
        if "sample" in report_df.columns:
            treatment_col = "sample"
        elif "Sample" in report_df.columns:
            treatment_col = "Sample"
        else:
            print(f"compare_treatments: no '{treatment_col}' or 'sample' column found. "
                  f"Pass a treatment_map to batch_report to add treatment annotations.")
            return report_df

    ctrl = report_df[report_df[treatment_col] == control_label]
    if len(ctrl) == 0:
        available = report_df[treatment_col].unique().tolist()
        print(f"compare_treatments: control '{control_label}' not found in "
              f"'{treatment_col}'. Available: {available}")
        return report_df

    ctrl_row = ctrl.iloc[0]
    numeric_cols = report_df.select_dtypes(include=[np.number]).columns

    result = report_df.copy()
    for col in numeric_cols:
        if ctrl_row[col] != 0 and not np.isnan(ctrl_row[col]):
            result[f"{col}_RR_vs_{control_label}"] = (
                (result[col] - ctrl_row[col]) / abs(ctrl_row[col]) * 100
            )
        else:
            result[f"{col}_RR_vs_{control_label}"] = np.nan

    return result
