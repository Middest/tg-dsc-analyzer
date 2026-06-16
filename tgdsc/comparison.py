"""Dual-atmosphere TG-DSC comparison: N2 vs Air.

Decouples SOM thermal stability (N2, pyrolysis only) from oxidative resistance
(Air, pyrolysis + combustion). The key interpretative framework:

  N2 T50  → intrinsic chemical bond strength of SOM
  Air T50 → bond strength + antioxidant protection from biochar/minerals
  Δ(Air - N2) → contribution of antioxidant protection mechanisms

Reference: Filimonenko et al. 2025 (Carbon Research).
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
import warnings


def compare_atmospheres(
    samples_n2: dict[str, pd.DataFrame],
    samples_air: dict[str, pd.DataFrame],
    treatment_map: dict[str, dict] | None = None,
    T_range: tuple[float, float] = (190, 640),
    pools_scheme: str = "filimonenko2025",
) -> dict[str, pd.DataFrame]:
    """Comprehensive N2 vs Air atmosphere comparison.

    Parameters
    ----------
    samples_n2, samples_air : dict mapping sample name to DataFrame.
    treatment_map : optional {name: {col: value}} annotations.
    T_range : temperature range for analysis.
    pools_scheme : thermal pool scheme.

    Returns
    -------
    dict with keys:
      - stability_comparison: side-by-side T50, TI, ED, Ea
      - pool_comparison: pool proportions by atmosphere
      - kinetics_comparison: Ea by pool by atmosphere
      - enhancement: Air/N2 ratios and Air-N2 deltas
    """
    from .stability import batch_report
    from .pools import pools_summary
    from .kinetics import pool_kinetics

    # Stability reports
    report_n2 = batch_report(samples_n2, treatment_map, pools_scheme, T_range)
    report_n2["Atmosphere"] = "N2"
    report_air = batch_report(samples_air, treatment_map, pools_scheme, T_range)
    report_air["Atmosphere"] = "Air"

    # Combine
    stability = pd.concat([report_n2, report_air], ignore_index=True)

    # Pool comparison
    pool_n2 = pools_summary(samples_n2, treatment_map, pools_scheme, T_range)
    pool_n2["Atmosphere"] = "N2"
    pool_air = pools_summary(samples_air, treatment_map, pools_scheme, T_range)
    pool_air["Atmosphere"] = "Air"
    pool_cmp = pd.concat([pool_n2, pool_air], ignore_index=True)

    # Kinetics comparison
    kin_frames = []
    for atm, samples in [("N2", samples_n2), ("Air", samples_air)]:
        for name, df in samples.items():
            kin = pool_kinetics(df, pools_scheme)
            kin["sample"] = name
            kin["Atmosphere"] = atm
            if treatment_map and name in treatment_map:
                for k, v in treatment_map[name].items():
                    kin[k] = v
            kin_frames.append(kin)
    kin_cmp = pd.concat(kin_frames, ignore_index=True)

    # Enhancement metrics
    enhance_rows = []
    sample_col = "sample" if "sample" in report_n2.columns else "Sample"
    if sample_col not in report_n2.columns:
        # Try to infer from treatment_map
        sample_col = None

    n2_samples = report_n2[sample_col].unique() if sample_col else report_n2.index
    air_samples = report_air[sample_col].unique() if sample_col else report_air.index

    if sample_col:
        n2_idx = report_n2.set_index(sample_col)
        air_idx = report_air.set_index(sample_col)
        common_samples = [s for s in n2_samples if s in air_samples]
    else:
        common_samples = list(range(len(report_n2)))

    key_metrics = ["T50_C", "ED_kJ_g_OM", "TI", "Ea_weighted_kJ_mol"]

    for s in common_samples:
        row = {"sample": s}
        if sample_col:
            n2_row = n2_idx.loc[s]
            air_row = air_idx.loc[s]
        else:
            n2_row = report_n2.iloc[s]
            air_row = report_air.iloc[s]

        row["SOM_loss_N2_pct"] = n2_row.get("SOM_loss_percent", np.nan)
        row["SOM_loss_Air_pct"] = air_row.get("SOM_loss_percent", np.nan)
        row["delta_SOM_loss_pct"] = row["SOM_loss_Air_pct"] - row["SOM_loss_N2_pct"]

        for m in key_metrics:
            n2_val = n2_row.get(m, np.nan)
            air_val = air_row.get(m, np.nan)
            row[f"{m}_N2"] = n2_val
            row[f"{m}_Air"] = air_val
            row[f"delta_{m}"] = air_val - n2_val if (
                not np.isnan(n2_val) and not np.isnan(air_val)
            ) else np.nan
            row[f"ratio_{m}"] = air_val / n2_val if (
                n2_val and not np.isnan(n2_val) and not np.isnan(air_val)
            ) else np.nan

        enhance_rows.append(row)

    enhancement = pd.DataFrame(enhance_rows)

    return {
        "stability_comparison": stability,
        "pool_comparison": pool_cmp,
        "kinetics_comparison": kin_cmp,
        "enhancement": enhancement,
    }


def atmosphere_shift_summary(
    comparison: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Extract the key Air-N2 shift metrics as a concise summary table."""
    enh = comparison["enhancement"]
    cols = ["sample", "delta_T50_C", "delta_TI", "delta_Ea_weighted_kJ_mol",
            "delta_ED_kJ_g_OM", "delta_SOM_loss_pct"]
    available = [c for c in cols if c in enh.columns]
    return enh[available].copy()


def atmosphere_shift_interpretation(
    delta_t50: float,
) -> dict:
    """Interpret the atmosphere T50 shift.

    Parameters
    ----------
    delta_t50 : Air T50 - N2 T50.

    Returns
    -------
    dict with interpretation and mechanim keywords.
    """
    if delta_t50 < -10:
        interpretation = (
            "Strongly negative shift: unprotected SOM decomposes at lower "
            "temperature in air due to O2-catalyzed oxidation. No significant "
            "antioxidant protection."
        )
        mechanism = "exposed_SOM_free_radical_oxidation"
    elif delta_t50 < -3:
        interpretation = (
            "Slightly negative shift: some SOM protection but O2 still "
            "lowers decomposition temperature."
        )
        mechanism = "weak_protection"
    elif delta_t50 < 5:
        interpretation = (
            "Minimal shift: SOM intrinsic thermal stability dominates; "
            "antioxidant protection roughly balances O2-catalysis."
        )
        mechanism = "balanced"
    elif delta_t50 < 20:
        interpretation = (
            "Moderate positive shift: biochar/mineral encapsulation provides "
            "measurable antioxidant shielding, requiring higher temperature "
            "for combustion."
        )
        mechanism = "organo_mineral_encapsulation"
    elif delta_t50 < 40:
        interpretation = (
            "Strong positive shift: significant antioxidant protection from "
            "biochar-derived condensed aromatics and organo-mineral complexes. "
            "O2 access to SOM is kinetically hindered."
        )
        mechanism = "biochar_antioxidant_shielding"
    else:
        interpretation = (
            "Very strong positive shift: extensive biochar encapsulation "
            "and/or mineral-biochar-SOM ternary complexes. SOM is highly "
            "resistant to oxidative thermal decomposition."
        )
        mechanism = "extensive_biochar_encapsulation"

    return {
        "delta_T50_C": delta_t50,
        "interpretation": interpretation,
        "mechanism": mechanism,
    }
