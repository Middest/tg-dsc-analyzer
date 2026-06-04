---
name: tg-dsc-analysis
description: >-
  TG-DSC thermal analysis for soil organic carbon stability. Covers TG/DTG curves, DSC heat flow,
  thermal pool partitioning (Filimonenko 2025), Coats-Redfern and iso-conversional kinetics,
  free DTG/DSC peak deconvolution, pool-constrained Gaussian fitting, and thermal stability
  indices (T50, ED, TI, R50). Use whenever the user asks about TG-DSC, thermogravimetric
  analysis, differential scanning calorimetry, thermal stability of SOM/biochar/soil, or
  needs to process NETZSCH STA data.
---

# TG-DSC Thermal Analysis Skill (v0.2.0)

Comprehensive TG-DSC data analysis toolkit for soil organic carbon thermal stability
characterization, powered by the `tg-dsc-analyzer` Python package.

## When to Use

- User mentions TG-DSC, TGA, DSC, thermogravimetric, differential scanning calorimetry
- User needs to analyze NETZSCH STA instrument data
- User wants thermal stability indices (T50, ED, TI, R50)
- User needs kinetics analysis (Coats-Redfern, KAS, FWO)
- User wants DTG/DSC peak deconvolution (free or pool-constrained)
- User is characterizing SOM, biochar, or soil thermal properties

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
from tgdsc.io import batch_read
from tgdsc.preprocess import smooth_savgol
from tgdsc.stability import batch_report, compare_treatments
from tgdsc.pools import pools_summary
from tgdsc.kinetics import pool_kinetics
from tgdsc.deconvolution import batch_deconvolve_pools
from tgdsc.visualize import *
from tgdsc.export import *

# 1. Load & preprocess
samples = batch_read(data_dir, pattern="ExpDat_*.xlsx", recursive=True)
samples = {k: smooth_savgol(v) for k, v in samples.items()}

# 2. Comprehensive thermal stability report
report = batch_report(samples, treatment_map, T_range=(190, 640))
report = compare_treatments(report, "CK")

# 3. Pool partitioning (Filimonenko 2025)
pool_df = pools_summary(samples, treatment_map)

# 4. Kinetics (Coats-Redfern per pool)
kin_df = pool_kinetics(samples)

# 5. Constrained pool deconvolution (v0.2.0)
deconv_df = batch_deconvolve_pools(samples)

# 6. Visualize
plot_tg_dtg(samples, T_range, save_path="output/tg_dtg.svg")
plot_pools_stacked(pool_df, save_path="output/pools.svg")
plot_constrained_deconvolution(fit_result, save_path="output/deconv.svg")
plot_pool_peak_comparison(deconv_df, metric="T_center_C", save_path="output/peaks.svg")

# 7. Export
to_excel({"Report": report, "Pools": pool_df, "Kinetics": kin_df, "Deconv": deconv_df},
         "output/full_report.xlsx")
```

### Phase 4: Deliverables

1. **Comprehensive report** — T50, ED, TI, SOM loss, weighted Ea per sample
2. **Pool analysis** — mass loss, proportion, ED per pool per sample
3. **Kinetics table** — Coats-Redfern Ea and R² per pool
4. **Constrained deconvolution** — peak center, width (FWHM), area per pool
5. **Figures** (SVG + PDF):
   - TG + DTG dual-panel
   - DSC curves
   - Pool stacked bar chart
   - Stability index comparison
   - Constrained deconvolution with pool shading
   - Pool peak center / FWHM comparison
6. **Export** — Excel (multi-sheet), CSV, JSON

## Deconvolution Methods

### Free Deconvolution (`deconvolve_dtg`)
- Fits N unconstrained peaks to DTG/DSC curve
- Use case: discovering unknown components
- Limitation: peak numbering inconsistent across samples

### Pool-Constrained Deconvolution (`deconvolve_pools`) — v0.2.0
- Fits exactly 1 Gaussian peak per thermal pool with center bounded to pool range
- Peak indices directly comparable: Peak 0 = labile, Peak 1 = stable, etc.
- Reveals intra-pool peak shifts (e.g., POC labile 210°C vs MAOC 244°C)
- Best used for: comparing peak center temperature and width between groups

## Thermal Pool Schemes

| Scheme | Labile | Stable | Persistent | Refractory |
|--------|--------|--------|------------|------------|
| filimonenko2025 | 190–390 | 390–490 | 490–590 | 590–640 |
| an_125_650 | 125–390 | 390–490 | 490–590 | 590–650 |

## Supported Formats

| Format | Columns | Origin |
|--------|---------|--------|
| NETZSCH Proteus | 6 col (Temp/Time/DSC/Mass/DTG/Sens) | STA 449F5 export |
| NETZSCH Legacy | 9 col (3x Time/Temp/Signal) | Older STA export |
| Generic CSV | Auto-detected | Any |

Unit conversions handled automatically: DSC mW/mg → μW/mg, DTG %/min → %/°C.

## Validation

Validated against manual Excel calculations (Haicheng paddy soil; CK/BC7.5/BC15/BC30):
- T50: ±0.02%
- Mass loss: ±0.08%
- ED: ±0.03%
- Ea: ±1.3%

PG-POC MAOC dataset (43 samples, 2026): constrained deconvolution reveals +34°C labile peak shift in MAOC vs POC.

## Package Location

```
C:/Users/Administrator/Desktop/tg-dsc-analyzer/
```

GitHub: https://github.com/Middest/tg-dsc-analyzer
