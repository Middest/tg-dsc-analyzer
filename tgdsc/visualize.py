"""Visualization: publication-quality TG-DSC figures.

Default style targets Nature-family journal requirements:
  - Clean sans-serif fonts
  - 7-8 pt font sizes
  - Multi-panel layout support
  - SVG + PDF export
"""

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


# --- global style ---

def set_style():
    """Apply Nature-style defaults."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 6.5,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.minor.width": 0.3,
        "ytick.minor.width": 0.3,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.minor.size": 2,
        "ytick.minor.size": 2,
        "lines.linewidth": 0.8,
        "figure.dpi": 300,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    })


COLORS = {
    "CK": "#2c3e50",
    "BC7.5": "#3498db",
    "BC15": "#e67e22",
    "BC30": "#e74c3c",
}

POOL_FILLS = {
    "labile": "#fdebd0",
    "stable": "#d5f5e3",
    "persistent": "#d6eaf8",
    "refractory": "#e8daef",
}


def _save(fig, save_path: str | Path | None, formats: list[str] | None = None):
    """Save figure to specified formats."""
    if save_path is None:
        return
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    if formats is None:
        formats = ["svg", "pdf"]
    for fmt in formats:
        fpath = save_path.with_suffix(f".{fmt}")
        fig.savefig(fpath, format=fmt)
    plt.close(fig)


# --- plot functions ---

def plot_tg(
    samples: dict[str, pd.DataFrame],
    treatment_map: dict[str, dict] | None = None,
    T_range_highlight: tuple[float, float] | None = (190, 640),
    save_path: str | Path | None = None,
    figsize: tuple[float, float] = (4, 3),
) -> plt.Figure:
    """TG (mass loss) curves for all samples."""
    set_style()
    fig, ax = plt.subplots(figsize=figsize)

    for name, df in samples.items():
        label = name
        color = COLORS.get(name, None)
        ax.plot(df["Temp_C"], df["TG_percent"], label=label, color=color, lw=0.8)

    if T_range_highlight:
        ax.axvspan(T_range_highlight[0], T_range_highlight[1], alpha=0.08, color="gray")

    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("TG (%)")
    ax.legend(frameon=False, loc="lower left")
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    _save(fig, save_path)
    return fig


def plot_dtg(
    samples: dict[str, pd.DataFrame],
    T_range: tuple[float, float] = (190, 640),
    save_path: str | Path | None = None,
    figsize: tuple[float, float] = (4, 3),
    smoothing_window: int | None = 15,
) -> plt.Figure:
    """DTG (derivative thermogravimetry) curves."""
    set_style()
    fig, ax = plt.subplots(figsize=figsize)

    for name, df in samples.items():
        color = COLORS.get(name, None)
        dtg = df["DTG_percent_per_C"].values
        T = df["Temp_C"].values
        if smoothing_window:
            from scipy.signal import savgol_filter
            dtg = savgol_filter(dtg, smoothing_window, 3)
        ax.plot(T, dtg, label=name, color=color, lw=0.8)

    ax.set_xlim(T_range)
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("DTG (%/°C)")
    ax.legend(frameon=False)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    _save(fig, save_path)
    return fig


def plot_tg_dtg(
    samples: dict[str, pd.DataFrame],
    T_range: tuple[float, float] = (190, 640),
    save_path: str | Path | None = None,
    figsize: tuple[float, float] = (6, 4),
) -> plt.Figure:
    """Combined TG + DTG dual-panel figure."""
    set_style()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)

    for name, df in samples.items():
        color = COLORS.get(name, None)
        ax1.plot(df["Temp_C"], df["TG_percent"], label=name, color=color, lw=0.8)

        from scipy.signal import savgol_filter
        dtg = savgol_filter(df["DTG_percent_per_C"].values, 15, 3)
        T = df["Temp_C"].values
        ax2.plot(T, dtg, color=color, lw=0.8)

    ax1.set_ylabel("TG (%)")
    ax1.legend(frameon=False, loc="lower left", ncol=2)
    ax1.set_xlim(T_range)
    ax1.axhline(100, ls="--", color="gray", lw=0.3, alpha=0.5)

    ax2.set_xlabel("Temperature (°C)")
    ax2.set_ylabel("DTG (%/°C)")

    for ax in [ax1, ax2]:
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    fig.tight_layout()
    _save(fig, save_path)
    return fig


def plot_dsc(
    samples: dict[str, pd.DataFrame],
    T_range: tuple[float, float] = (190, 640),
    save_path: str | Path | None = None,
    figsize: tuple[float, float] = (4, 3),
) -> plt.Figure:
    """DSC heat flow curves."""
    set_style()
    fig, ax = plt.subplots(figsize=figsize)

    for name, df in samples.items():
        color = COLORS.get(name, None)
        ax.plot(df["Temp_C"], df["DSC_uW_per_mg"], label=name, color=color, lw=0.8)

    ax.axhline(0, ls="--", color="gray", lw=0.3)
    ax.set_xlim(T_range)
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("DSC (μW/mg)")
    ax.legend(frameon=False)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    _save(fig, save_path)
    return fig


def plot_pools_stacked(
    pool_summary_df: pd.DataFrame,
    save_path: str | Path | None = None,
    figsize: tuple[float, float] = (5, 4),
    bar_width: float = 0.6,
) -> plt.Figure:
    """Stacked bar chart of thermal pool contributions."""
    set_style()
    pool_order = ["labile", "stable", "persistent", "refractory"]
    pool_labels = {
        "labile": "Labile (190–390°C)",
        "stable": "Stable (390–490°C)",
        "persistent": "Persistent (490–590°C)",
        "refractory": "Refractory (590–640°C)",
    }

    # Pivot to wide format
    samples = pool_summary_df["sample"].unique()
    pivot = pool_summary_df.pivot_table(
        index="sample", columns="pool", values="pool_proportion_percent",
        aggfunc="first",
    )

    # Reorder
    avail_pools = [p for p in pool_order if p in pivot.columns]
    pivot = pivot[avail_pools]

    fig, ax = plt.subplots(figsize=figsize)
    bottom = np.zeros(len(pivot))
    colors_list = [POOL_FILLS.get(p, "#cccccc") for p in avail_pools]

    for i, pool in enumerate(avail_pools):
        vals = pivot[pool].values
        ax.bar(range(len(pivot)), vals, bar_width, bottom=bottom,
               label=pool_labels.get(pool, pool), color=colors_list[i],
               edgecolor="white", lw=0.3)
        bottom += vals

    ax.set_xticks(range(len(pivot)))
    ax.set_xticklabels(pivot.index, rotation=0)
    ax.set_ylabel("Pool proportion (% of SOM)")
    ax.legend(frameon=False, fontsize=6)
    ax.set_ylim(0, 105)

    _save(fig, save_path)
    return fig


def plot_ea_comparison(
    kinetics_summary: pd.DataFrame,
    save_path: str | Path | None = None,
    figsize: tuple[float, float] = (5, 3.5),
) -> plt.Figure:
    """Grouped bar chart of Ea per pool per treatment."""
    set_style()
    fig, ax = plt.subplots(figsize=figsize)

    pool_cols = [c for c in kinetics_summary.columns if c.startswith("Ea_") and c.endswith("_kJ_mol")]
    samples = kinetics_summary.get("sample", kinetics_summary.index)

    x = np.arange(len(pool_cols))
    width = 0.8 / len(samples)

    for i, (_, row) in enumerate(kinetics_summary.iterrows()):
        vals = [row.get(c, np.nan) for c in pool_cols]
        offset = (i - len(samples) / 2 + 0.5) * width
        ax.bar(x + offset, vals, width, label=row.get("sample", i))

    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("Ea_", "").replace("_kJ_mol", "") for c in pool_cols])
    ax.set_ylabel("Ea (kJ/mol)")
    ax.legend(frameon=False)

    _save(fig, save_path)
    return fig


def plot_deconvolution(
    fit_result: dict,
    save_path: str | Path | None = None,
    figsize: tuple[float, float] = (5, 3.5),
    title: str = "",
) -> plt.Figure:
    """Show deconvolution fit with individual components."""
    set_style()
    fig, ax = plt.subplots(figsize=figsize)

    x = fit_result["x"]
    y = fit_result["y"]
    best_fit = fit_result["best_fit"]
    components = fit_result.get("components", [])
    result = fit_result.get("result")

    ax.plot(x, y, "k.", ms=1, alpha=0.3, label="Data")
    ax.plot(x, best_fit, "r-", lw=1, label=f"Fit (R²={fit_result.get('r2', 0):.3f})")

    # Plot each component if we have the result object
    if result is not None and components:
        colors = plt.cm.tab10(np.linspace(0, 1, len(components)))
        for i, comp in enumerate(components):
            prefix = f"p{i}_"
            comp_y = result.eval_components(x=x).get(prefix, None)
            if comp_y is not None:
                ax.plot(x, comp_y, "--", color=colors[i], lw=0.6,
                        label=f"Peak {i+1} ({comp['center']:.0f}°C)")

    if title:
        ax.set_title(title)
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Signal")
    ax.legend(frameon=False, fontsize=6)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))

    _save(fig, save_path)
    return fig


def plot_constrained_deconvolution(
    fit_result: dict,
    pools_scheme: str = "filimonenko2025",
    save_path: str | Path | None = None,
    figsize: tuple[float, float] = (6, 4.5),
    title: str = "",
) -> plt.Figure:
    """Plot constrained pool deconvolution with pool regions shaded."""
    set_style()
    from tgdsc.pools import get_pool_scheme

    pools = get_pool_scheme(pools_scheme)
    x = fit_result.get("x")
    y = fit_result.get("y")
    best_fit = fit_result.get("best_fit")
    components_y = fit_result.get("components_y", [])
    peaks = fit_result.get("peaks", [])

    if x is None or y is None:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        return fig

    fig, ax = plt.subplots(figsize=figsize)
    pool_colors_bg = ["#fdebd0", "#d5f5e3", "#d6eaf8", "#e8daef"]
    for i, p in enumerate(pools):
        ax.axvspan(p["T_low"], p["T_high"], alpha=0.25,
                   color=pool_colors_bg[i % len(pool_colors_bg)],
                   label="{} ({}-{}C)".format(p["pool"], p["T_low"], p["T_high"]))

    ax.plot(x, y, "k.", ms=1, alpha=0.3)
    comp_colors = ["#e67e22", "#27ae60", "#2980b9", "#8e44ad"]
    for i, comp_y in enumerate(components_y):
        ax.plot(x, comp_y, "--", color=comp_colors[i % len(comp_colors)], lw=1.0,
                label="Peak {} ({}C)".format(i+1, peaks[i]["T_center_C"]) if i < len(peaks) else None)

    if best_fit is not None:
        ax.plot(x, best_fit, "r-", lw=1.2, label="Fit (R²={:.3f})".format(fit_result.get("r2", 0)))

    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("DTG (%/°C)")
    if title:
        ax.set_title(title, fontsize=8)
    ax.legend(frameon=False, fontsize=5.5, loc="upper right")
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    _save(fig, save_path)
    return fig


def plot_pool_peak_comparison(
    deconv_df: pd.DataFrame,
    metric: str = "T_center_C",
    save_path: str | Path | None = None,
    figsize: tuple[float, float] = (8, 5),
) -> plt.Figure:
    """Compare constrained peak parameters across samples grouped by pool."""
    set_style()
    pools = deconv_df["pool"].unique()
    n_pools = len(pools)
    metric_labels = {
        "T_center_C": "Peak center (°C)",
        "fwhm_C": "Peak FWHM (°C)",
        "area": "Peak area",
        "area_fraction": "Area fraction",
        "r2": "Per-pool R²",
    }

    fig, axes = plt.subplots(1, n_pools, figsize=figsize, sharey=False)
    if n_pools == 1:
        axes = [axes]

    for i, pool in enumerate(pools):
        ax = axes[i]
        pool_data = deconv_df[deconv_df["pool"] == pool]
        vals = pool_data[metric].values
        names = pool_data["sample"].values
        colors_list = [COLORS.get(s, "#888888") for s in names]
        ax.bar(range(len(vals)), vals, color=colors_list, edgecolor="white", lw=0.3)
        ax.set_xticks(range(len(vals)))
        ax.set_xticklabels(names, rotation=90, ha="center", fontsize=5)
        ax.set_title(pool, fontsize=8, color=POOL_FILLS.get(pool, "#333"))
        ax.set_ylabel(metric_labels.get(metric, metric), fontsize=7)
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    fig.tight_layout()
    _save(fig, save_path)
    return fig


def plot_stability_comparison(
    report_df: pd.DataFrame,
    save_path: str | Path | None = None,
    figsize: tuple[float, float] = (7, 6),
) -> plt.Figure:
    """Multi-panel comparison of key stability indices."""
    set_style()
    metrics = [
        ("T50_C", "T50 (°C)"),
        ("ED_kJ_g_OM", "ED (kJ/g OM)"),
        ("SOM_loss_percent", "SOM loss (%)"),
        ("TI", "Thermostability Index"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    samples = report_df.get("sample", report_df.index)

    for (col, ylabel), ax in zip(metrics, axes.flat):
        if col not in report_df.columns:
            continue
        vals = report_df[col].values
        colors_list = [COLORS.get(s, "#888888") for s in samples]
        bars = ax.bar(range(len(vals)), vals, color=colors_list, edgecolor="white", lw=0.3)
        ax.set_xticks(range(len(vals)))
        ax.set_xticklabels(samples, rotation=30, ha="right", fontsize=6)
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel, fontsize=7)
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    fig.tight_layout()
    _save(fig, save_path)
    return fig
