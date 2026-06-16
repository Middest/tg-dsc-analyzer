"""Data I/O for TG-DSC instrument exports.

Supports two NETZSCH export formats:
  1. Legacy 9-column: Time|Temp|TG|Time|Temp|DTG|Time|Temp|DSC
  2. Proteus 6-column: Temp|Time|DSC(mW/mg)|Mass(%)|DTG(%/min)|Sensitivity

Both are normalized to: Time_min, Temp_C, TG_percent, DTG_percent_per_C, DSC_uW_per_mg.
"""

from pathlib import Path
import warnings
import numpy as np
import pandas as pd


def read_netzsch_excel(path: str | Path) -> pd.DataFrame:
    """Read a single NETZSCH STA Excel export (auto-detects format).

    Returns DataFrame with: Time_min, Temp_C, TG_percent,
    DTG_percent_per_C, DSC_uW_per_mg.
    """
    wb = pd.ExcelFile(path)
    ws = wb.parse(wb.sheet_names[0], header=None)

    ncols = ws.shape[1]

    if ncols == 9:
        return _read_legacy_9col(ws)
    elif ncols == 6:
        return _read_proteus(ws, path)
    else:
        return _read_generic(ws)


def _read_legacy_9col(ws: pd.DataFrame) -> pd.DataFrame:
    """Parse legacy 9-column NETZSCH export."""
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
    return out.dropna(subset=["Temp_C"]).reset_index(drop=True)


def _read_proteus(ws: pd.DataFrame, path: str | Path) -> pd.DataFrame:
    """Parse NETZSCH Proteus export (6 columns with metadata header).

    Format:
      Rows 0..N: metadata lines starting with '#'
      One header row starting with '##Temp.'
      Data below.

    Columns: Temp/°C, Time/min, DSC/(mW/mg), Mass/%, DTG/(%/min), Sensit./(uV/mW)
    Conversion: DSC mW/mg → uW/mg (×1000), DTG %/min → %/°C (÷ heating_rate).
    """
    # Find the '##Temp' header row
    for i in range(len(ws)):
        val = ws.iloc[i, 0]
        if isinstance(val, str) and "##Temp" in val:
            header_idx = i
            break
    else:
        # Fallback: try generic
        return _read_generic(ws)

    # Parse metadata for sample identity and instrument
    metadata = {}
    for i in range(header_idx):
        val = str(ws.iloc[i, 0]) if pd.notna(ws.iloc[i, 0]) else ""
        if val.startswith("#IDENTITY:"):
            metadata["identity"] = val.split(":", 1)[1].strip()
        elif val.startswith("#INSTRUMENT:"):
            metadata["instrument"] = val.split(":", 1)[1].strip()
        elif val.startswith("#MTYPE:"):
            metadata["measurement_type"] = val.split(":", 1)[1].strip()
        elif val.startswith("#DATE/TIME:"):
            metadata["datetime"] = val.split(":", 1)[1].strip()

    # Parse header row
    raw_headers = ws.iloc[header_idx].tolist()
    # Map columns by keyword
    col_indices = {}
    for j, h in enumerate(raw_headers):
        h_str = str(h).lower()
        if "temp" in h_str:
            col_indices["Temp_C"] = j
        elif "time" in h_str:
            col_indices["Time_min"] = j
        elif "dsc" in h_str:
            col_indices["DSC_raw"] = j
        elif "dtg" in h_str:
            col_indices["DTG_raw"] = j
        elif "mass" in h_str or "tg" in h_str:
            col_indices["TG_percent"] = j

    # Extract data
    data = ws.iloc[header_idx + 1:].copy()
    data = data.reset_index(drop=True)

    out = pd.DataFrame()

    for std_name, col_idx in col_indices.items():
        vals = pd.to_numeric(data.iloc[:, col_idx], errors="coerce")
        out[std_name] = vals

    # Convert units
    if "DSC_raw" in out.columns:
        # Proteus exports DSC in mW/mg, standard is uW/mg
        out["DSC_uW_per_mg"] = out.pop("DSC_raw") * 1000

    if "DTG_raw" in out.columns:
        # Proteus exports DTG in %/min, need %/°C
        # %/°C = %/min / heating_rate(°C/min)
        if "Temp_C" in out.columns and "Time_min" in out.columns:
            hr = _estimate_heating_rate(out["Temp_C"].values, out["Time_min"].values)
        else:
            hr = 10.0  # default
        out["DTG_percent_per_C"] = out.pop("DTG_raw") / hr

    # Ensure all standard columns
    for col in ["Time_min", "Temp_C", "TG_percent", "DTG_percent_per_C", "DSC_uW_per_mg"]:
        if col not in out.columns:
            out[col] = np.nan

    # Reorder
    out = out[["Time_min", "Temp_C", "TG_percent", "DTG_percent_per_C", "DSC_uW_per_mg"]]

    # Attach metadata as attrs
    out.attrs["metadata"] = metadata
    out.attrs["source"] = str(path)

    return out.dropna(subset=["Temp_C"]).reset_index(drop=True)


def _read_generic(ws: pd.DataFrame) -> pd.DataFrame:
    """Generic parser: find header row and detect columns."""
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


def _estimate_heating_rate(temp: np.ndarray, time: np.ndarray) -> float:
    """Estimate heating rate (C/min) from temperature-time data."""
    valid = ~np.isnan(temp) & ~np.isnan(time)
    if valid.sum() < 2:
        return 10.0
    t = time[valid]
    T = temp[valid]
    # Linear fit over 50-700 C range
    mask = (T > 50) & (T < 700)
    if mask.sum() > 10:
        slope, _ = np.polyfit(t[mask], T[mask], 1)
        return abs(slope)
    return abs((T[-1] - T[0]) / (t[-1] - t[0]))


def read_csv(path: str | Path) -> pd.DataFrame:
    """Read a generic CSV with TG-DSC columns."""
    df = pd.read_csv(path)
    col_map = _detect_columns(df.columns.tolist())
    return df.rename(columns=col_map)


def batch_read(
    folder: str | Path,
    pattern: str = "ExpDat_*.xlsx",
    reader: str = "netzsch_excel",
    file_filter: callable = None,
    recursive: bool = True,
) -> dict[str, pd.DataFrame]:
    """Batch-read all matching files in a folder.

    Parameters
    ----------
    folder : path to search.
    pattern : glob pattern for files.
    reader : "netzsch_excel" or "csv".
    file_filter : optional callable to filter filenames.
    recursive : if True, search subdirectories.

    Returns
    -------
    dict mapping sample name to DataFrame.
    """
    folder = Path(folder)
    if recursive:
        files = sorted(folder.rglob(pattern))
    else:
        files = sorted(folder.glob(pattern))
    if file_filter:
        files = [f for f in files if file_filter(f.name)]

    read_fn = {"netzsch_excel": read_netzsch_excel, "csv": read_csv}[reader]

    results = {}
    seen_names = set()
    for f in files:
        # Use filename stem (strip "ExpDat_" prefix) as sample ID
        if f.name.startswith("ExpDat_"):
            name = f.stem.replace("ExpDat_", "")
        else:
            name = f.stem

        # If all files are in the same directory, ensure unique names
        if name in seen_names:
            # Try parent_dir/filename to disambiguate
            name = f"{f.parent.name}_{name}"
        if name in seen_names:
            # Last resort: full path hash
            name = f"{name}_{abs(hash(str(f))) % 10000}"
        seen_names.add(name)

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
    """Merge multiple sample DataFrames into a long-format table."""
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
                for orig_key, orig_col in lower_map.items():
                    if all(k in orig_key for k in cand):
                        mapping[orig_col] = target
                        break
            elif cand in lower_map:
                mapping[lower_map[cand]] = target
                break

    return mapping
