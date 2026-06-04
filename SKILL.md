---
name: tg-dsc-analysis
description: >-
  TG-DSC thermal analysis for soil organic carbon stability. Covers TG/DTG curves, DSC heat flow,
  thermal pool partitioning (Filimonenko 2025), Coats-Redfern and iso-conversional kinetics,
  DTG/DSC peak deconvolution, and thermal stability indices (T50, ED, TI, R50).
  Use whenever the user asks about TG-DSC, thermogravimetric analysis, differential scanning
  calorimetry, thermal stability of SOM/biochar/soil, or needs to process NETZSCH STA data.
---

# TG-DSC Thermal Analysis Skill

Comprehensive TG-DSC data analysis toolkit for soil organic carbon thermal stability
characterization, powered by the `tg-dsc-analyzer` Python package.

## When to Use

- User mentions TG-DSC, TGA, DSC, thermogravimetric, differential scanning calorimetry
- User needs to analyze NETZSCH STA instrument data
- User wants thermal stability indices (T50, ED, TI, R50)
- User needs kinetics analysis (Coats-Redfern, KAS, FWO)
- User wants DTG/DSC peak deconvolution
- User is characterizing SOM, biochar, or soil thermal properties

## Setup (first use)

```bash
pip install -e "C:/Users/Administrator/Desktop/tg-dsc-analyzer"
```

## Workflow

### Phase 1: Data Discovery

Ask the user for the data directory path. Scan for:
- NETZSCH Excel exports (`.xlsx` with 9-column Time/Temp/TG/DTG/DSC format)
- Generic CSV files with temperature and signal columns
- `.sta` raw files (limited support)

Report back: number of samples found, temperature range, sample identifiers.

### Phase 2: Configuration

Confirm analysis parameters with the user:

| Parameter | Default | Description |
|-----------|---------|-------------|
| T_range | (190, 640) | SOM combustion range (Filimonenko 2025) |
| Pools scheme | filimonenko2025 | Labile/Stable/Persistent/Refractory boundaries |
| Kinetics method | Coats-Redfern | Single-rate or iso-conversional (needs multi-rate data) |
| Deconvolution | off | Whether to fit DTG/DSC peaks |
| Treatment mapping | user-provided | {sample_name: {Treatment, Biochar_rate, ...}} |

### Phase 3: Analysis

Run the pipeline:

```python
from tgdsc.io import batch_read, merge_samples
from tgdsc.preprocess import smooth_savgol
from tgdsc.stability import batch_report, compare_treatments
from tgdsc.pools import pools_summary
from tgdsc.kinetics import pool_kinetics, weighted_ea
from tgdsc.visualize import *
from tgdsc.export import *

# 1. Load
samples = batch_read(data_dir, pattern="*.xlsx", reader="netzsch_excel")
samples = {k: smooth_savgol(v) for k, v in samples.items()}

# 2. Analyze
report = batch_report(samples, treatment_map, pools_scheme, T_range)
report = compare_treatments(report, "CK")

# 3. Visualize & Export
plot_tg_dtg(samples, T_range, save_path="output/tg_dtg.svg")
plot_pools_stacked(pools_summary(samples, treatment_map), save_path="output/pools.svg")
to_excel({"Report": report, "Merged": merged}, "output/full_report.xlsx")
```

### Phase 4: Deliverables

The skill produces:
1. **Comprehensive report table** — T50, ED, TI, SOM loss, Ea per pool
2. **Pool analysis table** — mass loss, proportion, ED per thermal pool
3. **Kinetics table** — Ea and R² per pool per sample
4. **Figures** (SVG + PDF):
   - TG + DTG dual-panel plot
   - DSC heat flow curves
   - Pool stacked bar chart
   - Ea comparison bar chart
   - Stability index comparison
   - DTG/DSC deconvolution (if enabled)
5. **Export files** — Excel (.xlsx), CSV (.csv), JSON (.json)

## Key Methods Reference

### Thermal Pool Schemes

| Scheme | Labile | Stable | Persistent | Refractory |
|--------|--------|--------|------------|------------|
| filimonenko2025 | 190-390 | 390-490 | 490-590 | 590-640 |
| an_125_650 | 125-390 | 390-490 | 490-590 | 590-650 |

### Stability Indices

| Index | Formula | Meaning |
|-------|---------|---------|
| T50 | T at 50% cum. mass loss | Higher = more stable |
| ED | DSC energy / TG mass loss (kJ/g OM) | Energy stored per unit OM |
| TI | ML(350-550)/ML(200-550) | Higher = more recalcitrant |
| Ea | Coats-Redfern activation energy | Decomposition energy barrier |

### Kinetics Methods

- **Coats-Redfern (single-rate)**: First-order reaction, ln(-ln(1-α)/T²) vs 1/T
- **KAS (multi-rate)**: ln(β/T²) = const - Ea/(RT), needs ≥3 heating rates
- **FWO (multi-rate)**: ln(β) = const - 1.052·Ea/(RT), needs ≥3 heating rates

## Validation

For NETZSCH STA data from Haicheng soil study (CK/BC7.5/BC15/BC30):
- T50: ±0.02% vs Excel manual calculation
- Mass loss: ±0.08% vs Excel
- ED: ±0.03% vs Excel
- Ea: ±1.3% vs Excel

## Package Location

```
C:/Users/Administrator/Desktop/tg-dsc-analyzer/
```

GitHub: https://github.com/Middest/tg-dsc-analyzer
