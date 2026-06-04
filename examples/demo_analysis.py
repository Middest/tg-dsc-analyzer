"""Complete TG-DSC analysis demo using Haicheng soil data.

This script demonstrates the full analysis pipeline:
1. Data loading
2. TG/DTG/DSC analysis
3. Thermal pool characterization
4. Kinetics (Coats-Redfern)
5. Thermal stability indices
6. Visualization
7. Export
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from tgdsc.io import read_netzsch_excel, merge_samples
from tgdsc.preprocess import smooth_savgol
from tgdsc.tg import T50, mass_loss
from tgdsc.dsc import energy_density, detect_peaks
from tgdsc.pools import calculate_pools, pools_summary
from tgdsc.kinetics import pool_kinetics, weighted_ea
from tgdsc.stability import (
    thermostability_index, comprehensive_report, batch_report, compare_treatments,
)
from tgdsc.visualize import (
    plot_tg, plot_dtg, plot_dsc, plot_tg_dtg,
    plot_pools_stacked, plot_ea_comparison, plot_deconvolution,
    plot_stability_comparison, set_style,
)
from tgdsc.export import to_excel, summary_table, to_json

# --- Configuration ---
DATA_DIR = r"C:\Users\Administrator\Desktop\试验数据\2025\海城-2025\2025海城土壤\TG-DSC-海城"
OUTPUT_DIR = r"C:\Users\Administrator\Desktop\tg-dsc-analyzer\examples\output"

# File mapping (sample_id -> filename)
FILE_MAP = {
    "CK": "3.xlsx",
    "BC7.5": "4.xlsx",
    "BC15": "5.xlsx",
    "BC30": None,  # will be auto-detected (6-重复.xlsx)
}

TREATMENT_MAP = {
    "CK": {"Treatment": "CK", "Biochar_rate": 0},
    "BC7.5": {"Treatment": "BC7.5", "Biochar_rate": 7.5},
    "BC15": {"Treatment": "BC15", "Biochar_rate": 15},
    "BC30": {"Treatment": "BC30", "Biochar_rate": 30},
}

T_RANGE = (190, 640)
POOLS_SCHEME = "filimonenko2025"


def find_bc30_file(data_dir: str) -> str:
    """Find BC30 file (6-*.xlsx, excluding 6.xlsx)."""
    for f in sorted(os.listdir(data_dir)):
        if f.startswith("6-") and f.endswith(".xlsx") and "merged" not in f.lower():
            return f
    raise FileNotFoundError("BC30 file not found")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    set_style()

    # 1. Load data
    print("Loading data...")
    samples = {}
    for sid, fname in FILE_MAP.items():
        if fname is None:
            fname = find_bc30_file(DATA_DIR)
        samples[sid] = read_netzsch_excel(os.path.join(DATA_DIR, fname))
    print(f"  Loaded {len(samples)} samples: {list(samples.keys())}")

    # 2. Smooth
    samples = {k: smooth_savgol(v) for k, v in samples.items()}

    # 3. Comprehensive report
    print("\nCalculating comprehensive thermal stability report...")
    report_df = batch_report(samples, TREATMENT_MAP, POOLS_SCHEME, T_RANGE)
    report_df = compare_treatments(report_df, "CK")
    print(report_df[[
        "sample", "Treatment", "T50_C", "ED_kJ_g_OM",
        "SOM_loss_percent", "TI", "Ea_weighted_kJ_mol",
    ]].to_string(index=False))

    # 4. Pool analysis
    print("\nPool analysis...")
    pool_df = pools_summary(samples, TREATMENT_MAP, POOLS_SCHEME, T_RANGE)
    print(pool_df[["sample", "pool", "mass_loss_percent",
                    "pool_proportion_percent", "ED_kJ_g_OM"]].to_string(index=False))

    # 5. Kinetics per pool
    print("\nKinetics per pool...")
    for sid, df in samples.items():
        kin = pool_kinetics(df, POOLS_SCHEME)
        print(f"\n  {sid}:")
        print(kin[["pool", "Ea_kJ_mol", "R2", "n_points"]].to_string(index=False))

    # 6. Visualizations
    print("\nGenerating figures...")
    merged = merge_samples(samples, TREATMENT_MAP)

    plot_tg_dtg(samples, T_RANGE, save_path=f"{OUTPUT_DIR}/tg_dtg_combined.svg")
    print("  Saved TG-DTG combined figure")

    plot_dsc(samples, T_RANGE, save_path=f"{OUTPUT_DIR}/dsc_curves.svg")
    print("  Saved DSC figure")

    plot_pools_stacked(pool_df, save_path=f"{OUTPUT_DIR}/pools_stacked.svg")
    print("  Saved pools stacked bar chart")

    plot_stability_comparison(report_df, save_path=f"{OUTPUT_DIR}/stability_comparison.svg")
    print("  Saved stability comparison figure")

    # 7. Deconvolution demo (CK sample DTG)
    from tgdsc.deconvolution import deconvolve_dtg
    dtg_fit = deconvolve_dtg(samples["CK"], T_range=T_RANGE, n_components=4)
    if "error" not in dtg_fit:
        plot_deconvolution(dtg_fit, save_path=f"{OUTPUT_DIR}/dtg_deconvolution.svg",
                           title="CK DTG Deconvolution")
        print("  Saved DTG deconvolution figure")

    # 8. Export
    print("\nExporting results...")
    paths = {}

    # Excel workbook with all sheets
    paths["excel"] = to_excel({
        "Comprehensive_Report": report_df,
        "Pool_Analysis": pool_df,
        "Merged_Data": merged,
    }, f"{OUTPUT_DIR}/tg_dsc_full_report.xlsx")
    print(f"  Excel: {paths['excel']}")

    # CSV summary
    paths["csv"] = summary_table(
        {r["sample"]: r.to_dict() for _, r in report_df.iterrows()},
        f"{OUTPUT_DIR}/tg_dsc_summary.csv",
    )
    print(f"  CSV: {paths['csv']}")

    print("\n=== Analysis complete ===")
    print(f"Output directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
