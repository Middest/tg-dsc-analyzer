# tg-dsc-analyzer

TG-DSC thermal analysis toolkit for soil organic carbon stability characterization.

## Features

- **Data I/O**: Read NETZSCH STA Excel exports, CSV, batch multi-file import
- **TG/DTG Analysis**: Mass loss curves, characteristic temperatures, fractional conversion
- **DSC Analysis**: Heat flow curves, peak detection, energy integration
- **Thermal Pools**: Filimonenko 2025 scheme (Labile/Stable/Persistent/Refractory), customizable
- **Kinetics**: Coats-Redfern (single-rate), KAS and FWO (multi-rate iso-conversional)
- **Peak Deconvolution**: Gaussian/Lorentzian mixture model fitting for DTG and DSC
- **Stability Indices**: T50, ED (energy density), TI (thermostability index), R50
- **Visualization**: Nature-journal-style figures (SVG + PDF export)
- **Export**: Multi-sheet Excel, CSV, JSON

## Installation

```bash
pip install -e .
# or
pip install git+https://github.com/Middest/tg-dsc-analyzer.git
```

## Quick Start

```python
from tgdsc.io import batch_read
from tgdsc.preprocess import smooth_savgol
from tgdsc.stability import batch_report, compare_treatments
from tgdsc.visualize import plot_tg_dtg, plot_pools_stacked
from tgdsc.export import to_excel

# Load NETZSCH Excel exports
samples = batch_read("path/to/data/", pattern="*.xlsx")

# Smooth
samples = {k: smooth_savgol(v) for k, v in samples.items()}

# Map treatments
treatment_map = {
    "sample1": {"Treatment": "CK", "Biochar_rate": 0},
    "sample2": {"Treatment": "BC", "Biochar_rate": 15},
}

# Run full analysis
report = batch_report(samples, treatment_map)
report = compare_treatments(report, "CK")
print(report[["sample", "T50_C", "ED_kJ_g_OM", "TI", "Ea_weighted_kJ_mol"]])

# Visualize
plot_tg_dtg(samples, (190, 640), save_path="output/tg_dtg.svg")

# Export
to_excel({"Report": report}, "output/results.xlsx")
```

## Validation

Validated against manual Excel calculations using Haicheng paddy soil TG-DSC data:

| Metric | Deviation vs Excel |
|--------|-------------------|
| T50 | < 0.02% |
| Mass loss | < 0.08% |
| ED | < 0.03% |
| TI | < 0.20% |
| Ea | < 1.3% |

## Requirements

- Python >= 3.9
- numpy, scipy, pandas, matplotlib, openpyxl, lmfit

## License

MIT

## Citation

If you use this tool in your research, please cite:

> Xia F. (2026). tg-dsc-analyzer: TG-DSC thermal analysis toolkit for soil organic carbon. https://github.com/Middest/tg-dsc-analyzer

Method reference:
> Filimonenko, E. et al. (2025). SOM combustion range 190-640 °C with four-pool thermal stability partitioning.
