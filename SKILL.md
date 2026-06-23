---
name: tg-dsc-analysis
description: >-
  TG-DSC thermal analysis for soil organic carbon stability. Covers TG/DTG curves, DSC heat flow,
  thermal pool partitioning (Filimonenko 2025), Coats-Redfern and iso-conversional kinetics,
  free DTG/DSC peak deconvolution, pool-constrained Gaussian fitting, N2-vs-Air atmosphere
  comparison, and thermal stability indices (T50, ED, TI, R50). Use whenever the user asks
  about TG-DSC, thermogravimetric analysis, differential scanning calorimetry, thermal
  stability of SOM/biochar/soil, or needs to process NETZSCH STA data.
---

# TG-DSC Thermal Analysis Skill (v0.4.0)

Comprehensive TG-DSC data analysis toolkit for soil organic carbon thermal stability
characterization, powered by the `tg-dsc-analyzer` Python package.

## When to Use

- User mentions TG-DSC, TGA, DSC, thermogravimetric, differential scanning calorimetry
- User needs to analyze NETZSCH STA instrument data
- User wants thermal stability indices (T50, ED, TI, R50)
- User needs kinetics analysis (Coats-Redfern, KAS, FWO)
- User wants DTG/DSC peak deconvolution (free or pool-constrained)
- User is characterizing SOM, biochar, or soil thermal properties
- User needs N2 vs Air dual-atmosphere comparison

## Setup

```bash
pip install -e "C:/Users/Administrator/Desktop/试验数据/2026/2026-6 数据汇总/tg-dsc-analyzer"
```

> **Note**: Use `python` (not `python3`) in this environment — `python3` maps to the
> Windows App execution alias and fails. The package is installed at the path above
> (v0.3.0, installed as editable).

## Workflow

### Phase 1: Data Discovery

**Haicheng paddy soil data locations** (as of 2026-06):

| Atmosphere | Path | File pattern |
|------------|------|-------------|
| **Air** | `C:/Users/Administrator/Desktop/TG-HC-Air/` | `HC{n}/ExpDat_HC{n}.xlsx` |
| **N2** | `C:/Users/Administrator/Desktop/试验数据/2026/其他数据/数据/海城-gpt相关/TG/` | `{n}-TG.xlsx` |

**Treatment mapping** (verified 2026-06-22):
- HC5 / 5 = **CK** (control)
- HC6 / 6 = **BC7.5** (7.5 t/ha biochar)
- HC7 / 7 = **BC15** (15 t/ha biochar)
- HC8 / 8 = **BC30** (30 t/ha biochar)

**Data hierarchy** (bulk → fractions):
- `HC5`–`HC8` = bulk soil (全土)
- `HCPOC-5`–`HCPOC-8` = particulate organic carbon (颗粒态有机碳) from corresponding bulk
- `HCMAOC-5`–`HCMAOC-8` = mineral-associated organic carbon (矿质结合态有机碳) from corresponding bulk
- `HC1`, `HC2`, `HCPOC-1/2`, `HCMAOC-1/2` = additional samples

**TG-HC-Air directory** contains all 18 samples (HC1–HC8, HCPOC-1–8, HCMAOC-1–8).

Report: sample count, temperature range, instrument metadata.

**IMPORTANT**: If all export files are in a single flat directory (not per-sample subfolders),
`batch_read` will use the filename stem (stripping "ExpDat_" prefix) as the sample name.
If duplicate names would occur, disambiguation suffixes are appended automatically.
For complete control, load files individually with `read_netzsch_excel(path)`.

### Phase 2: Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| T_range | (190, 640) | SOM combustion range (Filimonenko 2025) |
| Pools scheme | filimonenko2025 | Labile/Stable/Persistent/Refractory |
| Kinetics | Coats-Redfern | Single-rate; KAS/FWO need multi-rate data |
| Deconvolution | pool-constrained | Free (n-peak) or constrained (1 peak per pool) |
| Treatment mapping | auto-group | By filename prefix |

### Phase 3: Analysis Pipeline

```python
from tgdsc.io import batch_read, read_netzsch_excel
from tgdsc.preprocess import smooth_savgol
from tgdsc.stability import batch_report, compare_treatments
from tgdsc.pools import pools_summary
from tgdsc.kinetics import pool_kinetics
from tgdsc.deconvolution import batch_deconvolve_pools
from tgdsc.export import to_excel

# ── 1. Load & preprocess ─────────────────────────────────────
# N2 data: flat dir, files named "5-TG.xlsx" ~ "8-TG.xlsx"
# Air data: per-sample subdirs, files named "HC5/ExpDat_HC5.xlsx"
codes = {"5": "CK", "6": "BC7.5", "7": "BC15", "8": "BC30"}
n2_dir  = r"C:\Users\Administrator\Desktop\试验数据\2026\其他数据\数据\海城-gpt相关\TG"
air_dir = r"C:\Users\Administrator\Desktop\TG-HC-Air"

samples_n2  = {treat: smooth_savgol(read_netzsch_excel(f"{n2_dir}/{code}-TG.xlsx"))
               for code, treat in codes.items()}
samples_air = {treat: smooth_savgol(read_netzsch_excel(f"{air_dir}/HC{code}/ExpDat_HC{code}.xlsx"))
               for code, treat in codes.items()}

# ── 2. Treatment map — values MUST be dicts, not strings ────
treatment_map = {t: {"Treatment": t} for t in codes.values()}

# ── 3. Thermal stability reports ─────────────────────────────
T_RANGE = (190, 640)  # Filimonenko 2025 (STOTEN) scheme
report_air = batch_report(samples_air, treatment_map, T_range=T_RANGE)
report_n2  = batch_report(samples_n2,  treatment_map, T_range=T_RANGE)
report_air = compare_treatments(report_air, "CK")
report_n2  = compare_treatments(report_n2,  "CK")

# ── 4. Pool partitioning (Filimonenko 2025) ──────────────────
pool_air = pools_summary(samples_air, treatment_map)
pool_n2  = pools_summary(samples_n2,  treatment_map)

# ── 5. Kinetics — pool_kinetics() takes SINGLE DataFrame ─────
# Must iterate; it does NOT accept a dict of samples
kin_rows = []
for samples, atm in [(samples_n2, "N2"), (samples_air, "Air")]:
    for name, df in samples.items():
        k = pool_kinetics(df)
        k["sample"] = name
        k["Atmosphere"] = atm
        kin_rows.append(k)
kin_df = pd.concat(kin_rows, ignore_index=True)

# ── 6. Constrained pool deconvolution ────────────────────────
deconv_air = batch_deconvolve_pools(samples_air)
deconv_n2  = batch_deconvolve_pools(samples_n2)
# Check fit_quality: "good" > 0.8, "marginal" 0.5-0.8, "poor_negative_r2" < 0
# NOTE: Column is fwhm_C (lowercase), not FWHM_C

# ── 7. Manual atmosphere shift calculation ───────────────────
# compare_atmospheres() may not work with mismatched treatment keys;
# manual calculation is more reliable:
rows = []
for t in codes.values():
    n2 = report_n2[report_n2["sample"] == t].iloc[0]
    air = report_air[report_air["sample"] == t].iloc[0]
    rows.append({"Treatment": t,
        "N2_T50_C": n2["T50_C"], "Air_T50_C": air["T50_C"],
        "delta_T50_C": air["T50_C"] - n2["T50_C"],
        "N2_SOM_loss_pct": n2["SOM_loss_percent"],
        "Air_SOM_loss_pct": air["SOM_loss_percent"],
        "delta_SOM_loss_pct": air["SOM_loss_percent"] - n2["SOM_loss_percent"],
    })
shift_df = pd.DataFrame(rows)

# ── 8. Export ────────────────────────────────────────────────
to_excel({
    "Stability_Report": pd.concat([report_air, report_n2], ignore_index=True),
    "Atmosphere_Shift": shift_df,
    "Pool_Summary": pd.concat([pool_air, pool_n2], ignore_index=True),
    "Kinetics": kin_df,
    "Deconvolution": pd.concat([deconv_air, deconv_n2], ignore_index=True),
}, "output/dual_atmosphere_report.xlsx")
```

### Phase 4: Deliverables

1. **Comprehensive report** — T50, ED, TI, SOM loss, weighted Ea per sample
2. **Pool analysis** — mass loss, proportion, ED per pool per sample
3. **Kinetics table** — Coats-Redfern Ea and R² per pool
4. **Constrained deconvolution** — peak center, width (FWHM), area, fit_quality per pool
5. **Figures** (SVG + PDF):
   - TG + DTG dual-panel
   - DSC curves
   - Pool stacked bar chart
   - Stability index comparison
   - Constrained deconvolution with pool shading
   - Pool peak center / FWHM comparison
6. **Dual-atmosphere comparison** — N2 vs Air stability metrics, enhancement ratios
7. **Export** — Excel (multi-sheet), CSV, JSON

## Deconvolution Methods

### Free Deconvolution (`deconvolve_dtg`)
- Fits N unconstrained peaks to DTG/DSC curve
- Use case: discovering unknown components
- Limitation: peak numbering inconsistent across samples

### Pool-Constrained Deconvolution (`deconvolve_pools`) — v0.3.0
- Fits exactly 1 Gaussian peak per thermal pool with center bounded to pool range
- Peak indices directly comparable: Peak 0 = labile, Peak 1 = stable, etc.
- Reveals intra-pool peak shifts (e.g., POC labile 210°C vs MAOC 244°C)
- Best used for: comparing peak center temperature and width between groups
- **v0.3.0**: Added `fit_quality` column — check this before interpreting results.
  Values: "good" (R²>0.8), "acceptable" (0.5-0.8), "marginal" (0-0.5),
  "poor_negative_r2" (<0, fit failed), "flat_signal" (DTG amplitude < 1e-5)

## Thermal Pool Schemes

| Scheme | Labile | Stable | Persistent | Refractory |
|--------|--------|--------|------------|------------|
| filimonenko2025 | 190–390 | 390–490 | 490–590 | 590–640 |
| an_125_650 | 125–390 | 390–490 | 490–590 | 590–650 |

## Dual-Atmosphere Interpretation (v0.3.0)

The N2 vs Air comparison decouples two dimensions of SOM stability:

| Metric | What it measures | Atmosphere |
|--------|-----------------|------------|
| N2 T50 | Intrinsic chemical bond strength | Pyrolysis only |
| Air T50 | Bond strength + antioxidant protection | Pyrolysis + combustion |
| ΔT50 (Air–N2) | Contribution of biochar/mineral antioxidant shielding | — |

Interpretation thresholds for ΔT50 = Air_T50 - N2_T50:
- **< −10°C**: Unprotected SOM, O2-catalyzed decomposition
- **−3 to +5°C**: Balanced, intrinsic stability dominates
- **+5 to +40°C**: Biochar/mineral antioxidant shielding active
- **> +40°C**: Extensive encapsulation, strong oxidative resistance

## Supported Formats

| Format | Columns | Origin |
|--------|---------|--------|
| NETZSCH Proteus | 6 col (Temp/Time/DSC/Mass/DTG/Sens) | STA 449F5 export |
| NETZSCH Legacy | 9 col (3x Time/Temp/Signal) | Older STA export |
| Generic CSV | Auto-detected | Any |

Unit conversions handled automatically: DSC mW/mg → μW/mg, DTG %/min → %/°C.

## Troubleshooting

### Negative ED (Energy Density)
If `ED_kJ_g_OM` is negative or extreme:
1. The DSC baseline likely drifted (DSC signal is not zero-centered)
2. Check `ED_kJ_g_OM_raw` for the unclamped value
3. Consider applying `baseline_correct_dsc()` before analysis
4. The `ed_warnings` field will contain diagnostic codes

### Deconvolution R² < 0
- DTG signal amplitude is too low for reliable peak fitting
- Common with low-SOM samples (<3% mass loss)
- The constrained deconvolution `fit_quality` column will indicate failure
- Alternative: use free deconvolution to determine resolvable peak count first

### batch_read returns fewer samples than expected
- All files in same directory: `batch_read` uses filenames, not folder names
- Check for name collisions in the output dictionary keys
- Use `read_netzsch_excel()` individually for full control

### Windows console encoding issues
- Unicode characters (✓, °C symbol) may display garbled in Windows CMD
- Use `[OK]` instead of ✓ in print statements
- Figures are unaffected (SVG/PDF use proper encoding)

## Validation

Validated against manual Excel calculations (Haicheng paddy soil; CK/BC7.5/BC15/BC30):
- T50: ±0.02%
- Mass loss: ±0.08%
- ED: ±0.03%
- Ea: ±1.3%

PG-POC MAOC dataset (43 samples, 2026): constrained deconvolution reveals +34°C labile peak shift in MAOC vs POC.

### N2 vs Air Dual-Atmosphere Results (2026-06-22, verified)

Bulk soil HC5–HC8 (CK, BC7.5, BC15, BC30), Filimonenko 2025 scheme:

| Treatment | T50 N2 | T50 Air | **ΔT50** | ΔSOM | Interpretation |
|-----------|--------|---------|----------|------|----------------|
| CK | 412.6 | 394.8 | **−17.8** | +0.44 | Unprotected SOM, O₂-catalyzed decomposition |
| BC7.5 | 413.0 | 408.6 | **−4.4** | +1.00 | Balanced, intrinsic stability dominates |
| BC15 | 415.1 | 442.0 | **+26.8** | +2.20 | Biochar antioxidant shielding active |
| BC30 | 415.3 | 439.7 | **+24.3** | +2.01 | Biochar antioxidant shielding active |

**Key insight**: Under N2, all treatments have similar T50 (~413–415°C) — intrinsic chemical
bond strength is similar. Under Air, a clear biochar dose effect emerges: CK shows O₂-catalyzed
decomposition (T50 drops in Air), while BC15/BC30 show strong antioxidant protection (T50 rises).
This confirms the **Binding-Mode Filter hypothesis** — biochar provides oxidative shielding
that decouples intrinsic stability from oxidative stability.

**Thermal pool redistribution (Air)**:
- CK: labile 48% → persistent 15%
- BC30: labile 36% → persistent 31%
- Biochar shifts SOM from labile toward persistent pool, with BC15 showing the strongest effect.

**Filimonenko 2025 (STOTEN) paper confirmation**: Zotero ID 1217, DOI `10.1016/j.scitotenv.2025.179934`.
Pool scheme: Labile 190–390, Stable 390–490, Persistent 490–590, Refractory 590–640 °C —
identical to the `filimonenko2025` scheme in the package.

## Known Issues

1. **BC30 N2 raw data anomaly**: The N2 TG data file (`8-TG.xlsx`) for BC30 produces a
   negative ED (−8.7 kJ/g in labile pool) with the standard pipeline — DSC baseline drift
   suspected. The validated Air data is unaffected. Use Air results for BC30 ED; for N2,
   the user's validated reference values (from `BC_paddy_OC_HC_by_class.xlsx`) should be
   used instead. Root cause: likely DSC baseline drift during N2 run; recommend
   re-exporting from NETZSCH.

2. **Low-SOM DTG deconvolution**: Bulk soil samples have low SOM (2.8–5.0% mass loss),
   producing flat DTG signals. Pool-constrained Gaussian fitting consistently yields
   `fit_quality = poor_negative_r2` (R² < 0) for all pools. Deconvolution results
   (peak centers, FWHM) are **unreliable for bulk soil** — use thermal stability
   indices (T50, ED, TI) and pool mass-loss partitioning instead. Deconvolution is
   more useful for POC/MAOC fractions with higher SOM content.

3. **Data file copy confusion**: The folder `TG相关数据/26060836483/` contains **copies
   of the Air data**, NOT the N2 data. The actual N2 data is at `试验数据/2026/其他数据/
   数据/海城-gpt相关/TG/` with files named `5-TG.xlsx` ~ `8-TG.xlsx`. Always verify
   data identity by comparing TG/DSC means before running dual-atmosphere analysis.

4. **`treatment_map` API**: Values must be dicts (`{"Treatment": "CK"}`), not plain
   strings. Passing plain strings causes `'str' object has no attribute 'items'` errors
   in `batch_report`.

5. **`pool_kinetics()` API**: Takes a **single** DataFrame, not a dict of samples.
   Must iterate manually over samples.

## Package Location

```
C:/Users/Administrator/Desktop/试验数据/2026/2026-6 数据汇总/tg-dsc-analyzer/
```

## Analysis Output Directory

```
C:/Users/Administrator/Desktop/试验数据/2026/2026-6 数据汇总/HC-TG分析结果/
```

GitHub: https://github.com/Middest/tg-dsc-analyzer
