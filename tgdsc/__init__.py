"""tg-dsc-analyzer: TG-DSC thermal analysis toolkit for soil organic carbon.

Submodules:
  io          — data import (NETZSCH Excel/Proteus, CSV, batch)
  preprocess  — smoothing, baseline correction
  tg          — TG/DTG analysis
  dsc         — DSC analysis
  pools       — thermal pool partitioning (Filimonenko 2025)
  kinetics    — Coats-Redfern, KAS, FWO
  deconvolution — free + pool-constrained Gaussian fitting
  stability   — T50, ED, TI, comprehensive report
  comparison  — dual-atmosphere (N2 vs Air) comparison
  visualize   — publication-quality figures
  export      — Excel/CSV/JSON export
"""

__version__ = "0.3.0"
