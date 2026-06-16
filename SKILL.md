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

# TG-DSC Thermal Analysis Skill (v0.3.0)

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
pip install -e "C:/Users/Administrator/Desktop/tg-dsc-analyzer"
```

## Workflow

### Phase 1: Data Discovery

Scan the data directory for:
- NETZSCH Proteus exports (`.xlsx`, 6-column with `##Temp.` header) — auto-detected
- NETZSCH legacy exports (`.xlsx`, 9-column) — auto-detected
- Generic CSV files

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
from tgdsc.comparison import compare_atmospheres, atmosphere_shift_interpretation
from tgdsc.visualize import *
from tgdsc.export import *

# 1. Load & preprocess (choose one approach)
# Option A: batch load (best when files are in per-sample subfolders)
samples = batch_read(data_dir, pattern="ExpDat_*.xlsx", recursive=True)
# Option B: load individually (best for flat directories or custom naming)
samples = {f"HC{i}": read_netzsch_excel(f"ExpDat_HC{i}.xlsx") for i in range(5,9)}
samples = {k: smooth_savgol(v) for k, v in samples.items()}

# 2. Comprehensive thermal stability report
report = batch_report(samples, treatment_map, T_range=(190, 640))
# NOTE: compare_treatments requires treatment_map for Treatment column;
# without it, it falls back to sample column or skips gracefully
report = compare_treatments(report, "CK")

# 3. Pool partitioning (Filimonenko 2025)
pool_df = pools_summary(samples, treatment_map)

# 4. Kinetics (Coats-Redfern per pool)
kin_df = pool_kinetics(samples)

# 5. Constrained pool deconvolution (v0.2.0)
# Check fit_quality column: "good" > 0.8, "marginal" 0.5-0.8, "poor_negative_r2" < 0
# If fit_quality is "flat_signal" or "poor_negative_r2", deconv results are unreliable
deconv_df = batch_deconvolve_pools(samples)

# 6. Dual-atmosphere comparison (v0.3.0)
comparison = compare_atmospheres(samples_n2, samples_air, treatment_map)
shift_df = comparison["enhancement"]  # Air-N2 deltas and ratios
for _, row in shift_df.iterrows():
    interp = atmosphere_shift_interpretation(row["delta_T50_C"])
    print(f"{row['sample']}: {interp['interpretation']}")

# 7. Visualize
plot_tg_dtg(samples, T_range, save_path="output/tg_dtg.svg")
plot_pools_stacked(pool_df, save_path="output/pools.svg")
plot_constrained_deconvolution(fit_result, save_path="output/deconv.svg")
plot_pool_peak_comparison(deconv_df, metric="T_center_C", save_path="output/peaks.svg")

# 8. Export
to_excel({"Report": report, "Pools": pool_df, "Kinetics": kin_df, "Deconv": deconv_df,
          "Enhancement": shift_df},
         "output/full_report.xlsx")
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

N2-vs-Air comparison (Haicheng paddy, 2026-06): BC30 shows +35.8°C T50 shift (Air > N2),
indicating extensive biochar antioxidant shielding at high application rates.

## Known Issues

1. **BC30 N2 raw data anomaly**: The NETZSCH export file `6.xlsx` / `8-TG.xlsx` for BC30
   under N2 produces a negative ED (−4.3 kJ/g) with the standard pipeline. The user's
   validated reference values (from `BC_paddy_OC_HC_by_class.xlsx`) should be used instead.
   Root cause likely DSC baseline drift; recommend re-exporting from NETZSCH.

2. **Low-SOM DTG deconvolution**: Samples with total SOM loss < 3% in the 190–640°C range
   produce very flat DTG signals; pool-constrained Gaussian fitting may fail (R² < 0).
   The `fit_quality` field will indicate this.

## Package Location

```
C:/Users/Administrator/Desktop/tg-dsc-analyzer/
```

GitHub: https://github.com/Middest/tg-dsc-analyzer
