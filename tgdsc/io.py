"""Data I/O for TG-DSC instrument exports.

Supports NETZSCH STA Excel exports, generic CSV, and multi-file batch import.
"""

from pathlib import Path
import warnings
import numpy as np
import pandas as pd


def read_netzsch_excel(path: str | Path) -> pd.DataFrame:
    """Read a single NETZSCH STA Excel export file.

    NETZSCH exports typically have 9 columns in 3 groups:
      Time | Temp | TG | Time | Temp | DTG | Time | Temp | DSC

    Returns a DataFrame with columns: Time_min, Temp_C, TG_percent,
    DTG_percent_per_C, DSC_uW_per_mg.
    """
    wb = pd.ExcelFile(path)
    ws = wb.parse(wb.sheet_names[0], header=None)

    # NETZSCH exports have 2 header rows then data
    if ws.shape[1] == 9:
        data_rows = ws.iloc[2:].copy()
        data_rows.columns = [
            "Time_min", "Temp_C", "TG_percent",
            "Time_min_dup", "Temp_C_dup", "DTG_percent_per_C",
            "Time_min_dup2", "Temp_C_dup2", "DSC_uW_per_mg",
        ]
        data_rows = data_rows.astype(float)
        out = pd.DataFrame({
            "Time_min": data_rows["Time_min"].values,
            "Temp_C": data_rows["Temp_C"].values,
            "TG_percent": data_rows["TG_percent"].values,
            "DTG_percent_per_C": data_rows["DTG_percent_per_C"].values,
            "DSC_uW_per_mg": data_rows["DSC_uW_per_mg"].values,
        })
    else:
        # Try single header row with fewer columns
        header_row = _find_header_row(ws)
        out = ws.iloc[header_row:].copy()
        out.columns = out.iloc[0]
        out = out.iloc[1:].reset_index(drop=True)
        col_map = _detect_columns(out.columns.tolist())
        out = out.rename(columns=col_map)
        for c in ["Time_min", "Temp_C", "TG_percent", "DTG_percent_per_C", "DSC_uW_per_mg"]:
            if c in out.columns:
                out[c] = pd.to_numeric(out[c], errors="coerce")

    return out.dropna(subset=["Temp_C"]).reset_index(drop=True)


def read_csv(path: str | Path) -> pd.DataFrame:
    """Read a generic CSV with TG-DSC columns."""
    df = pd.read_csv(path)
    col_map = _detect_columns(df.columns.tolist())
    return df.rename(columns=col_map)


def batch_read(
    folder: str | Path,
    pattern: str = "*.xlsx",
    reader: str = "netzsch_excel",
    file_filter: callable = None,
) -> dict[str, pd.DataFrame]:
    """Batch-read all matching files in a folder.

    Returns dict mapping sample name (stem) to DataFrame.
    """
    folder = Path(folder)
    files = sorted(folder.glob(pattern))
    if file_filter:
        files = [f for f in files if file_filter(f.name)]

    read_fn = {"netzsch_excel": read_netzsch_excel, "csv": read_csv}[reader]

    results = {}
    for f in files:
        name = f.stem
        try:
            df = read_fn(f)
            if len(df) > 0:
                results[name] = df
        except Exception as e:
            warnings.warn(f"Failed to read {f.name}: {e}")
    return results


def merge_samples(
    samples: dict[str, pd.DataFrame],
    treatment_map: dict[str, dict] | None = None,
) -> pd.DataFrame:
    """Merge multiple sample DataFrames into a long-format table.

    Parameters
    ----------
    samples : dict
        {sample_name: DataFrame} from batch_read.
    treatment_map : dict, optional
        {sample_name: {"Treatment": str, "Biochar_rate": float, ...}}

    Returns
    -------
    pd.DataFrame with columns: Sample_No, Treatment, Biochar_rate, ...,
    Time_min, Temp_C, TG_percent, DTG_percent_per_C, DSC_uW_per_mg.
    """
    frames = []
    for name, df in samples.items():
        sub = df.copy()
        sub["Sample_No"] = name
        sub["Source_file"] = name
        if treatment_map and name in treatment_map:
            for k, v in treatment_map[name].items():
                sub[k] = v
        frames.append(sub)

    merged = pd.concat(frames, ignore_index=True)
    col_order = [
        "Sample_No", "Treatment", "Biochar_rate", "Source_file",
        "Time_min", "Temp_C", "TG_percent", "DTG_percent_per_C", "DSC_uW_per_mg",
    ]
    existing = [c for c in col_order if c in merged.columns]
    extra = [c for c in merged.columns if c not in existing]
    return merged[existing + extra]


# --- helpers ---

def _find_header_row(ws: pd.DataFrame) -> int:
    """Find the first row that looks like a header."""
    for i, row in ws.iterrows():
        vals = [str(v).lower() for v in row if pd.notna(v)]
        if any(k in v for v in vals for k in ("temp", "time", "tg", "dtg", "dsc")):
            return i
    return 0


def _detect_columns(cols: list[str]) -> dict[str, str]:
    """Map fuzzy column names to standard names."""
    mapping = {}
    lower_map = {c.lower(): c for c in cols}

    patterns = {
        "time_min": [("time", "min"), "time_min", "time"],
        "temp_c": [("temp", "c"), "temp_c", "temperature", "temp"],
        "tg_percent": [("tg",), "tg_percent", "tg"],
        "dtg_percent_per_c": [("dtg",), "dtg_percent_per_c", "dtg"],
        "dsc_uw_per_mg": [("dsc",), "dsc_uw_per_mg", "dsc"],
    }

    for target, candidates in patterns.items():
        for cand in candidates:
            if isinstance(cand, tuple):
                # match columns that contain ALL keywords
                for orig_key, orig_col in lower_map.items():
                    if all(k in orig_key for k in cand):
                        mapping[orig_col] = target
                        break
            elif cand in lower_map:
                mapping[lower_map[cand]] = target
                break

    return mapping
