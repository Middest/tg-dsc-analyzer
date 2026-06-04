"""Export: data and results export to Excel, CSV, JSON."""

import json
from pathlib import Path
import numpy as np
import pandas as pd


def _sanitize(obj):
    """Convert numpy types to native Python for JSON serialization."""
    if isinstance(obj, (np.integer,)):
        return int(obj)
    elif isinstance(obj, (np.floating,)):
        return float(obj) if not np.isnan(obj) else None
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: _sanitize(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [_sanitize(v) for v in obj]
    return obj


def to_excel(
    data_dict: dict[str, pd.DataFrame],
    path: str | Path,
    merge_sheets: bool = True,
) -> Path:
    """Export multiple DataFrames to a multi-sheet Excel workbook.

    Parameters
    ----------
    data_dict : dict of {sheet_name: DataFrame}.
    path : output file path.
    merge_sheets : if True, when a single DataFrame is passed without a key,
                   still writes as a sheet.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if path.suffix not in (".xlsx",):
        path = path.with_suffix(".xlsx")

    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, df in data_dict.items():
            df.to_excel(writer, sheet_name=sheet_name[:31], index=False)  # Excel 31-char limit

    return path


def summary_table(
    comprehensive_reports: dict[str, dict],
    path: str | Path,
) -> Path:
    """Export comprehensive reports as a CSV summary."""
    path = Path(path)
    rows = []
    for sample, report in comprehensive_reports.items():
        row = {"sample": sample, **_sanitize(report)}
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    return path


def to_json(
    data: dict,
    path: str | Path,
) -> Path:
    """Export data to JSON."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix not in (".json",):
        path = path.with_suffix(".json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(_sanitize(data), f, indent=2, ensure_ascii=False)
    return path


def export_full_report(
    samples: dict[str, pd.DataFrame],
    pooled_samples: pd.DataFrame | None,
    report_df: pd.DataFrame,
    output_dir: str | Path,
    prefix: str = "tg_dsc",
) -> dict[str, Path]:
    """Export a complete analysis report package.

    Returns dict mapping description to file path.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths = {}

    # 1. Merged raw data
    if pooled_samples is not None:
        paths["merged_data"] = to_excel(
            {"Merged_for_stats": pooled_samples},
            output_dir / f"{prefix}_merged_data.xlsx",
        )

    # 2. Comprehensive report
    paths["report"] = summary_table(
        {r.get("sample", str(i)): r for i, (_, r) in enumerate(report_df.iterrows())},
        output_dir / f"{prefix}_comprehensive_report.csv",
    )

    # 3. Detailed Excel
    paths["excel"] = to_excel(
        {"Comprehensive_Report": report_df},
        output_dir / f"{prefix}_results.xlsx",
    )

    # 4. JSON
    paths["json"] = to_json(
        report_df.to_dict(orient="records"),
        output_dir / f"{prefix}_results.json",
    )

    return paths
