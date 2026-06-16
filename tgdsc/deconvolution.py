"""Peak deconvolution for DTG and DSC curves.

Fits overlapping peaks using Gaussian, Lorentzian, or Voigt mixture models
via lmfit. Essential for separating multi-component thermal degradation
profiles in complex organic materials (biochar-amended soils).
"""

import numpy as np
import pandas as pd
import warnings
from lmfit import Model
from lmfit.models import GaussianModel, LorentzianModel, VoigtModel


MODEL_MAP = {
    "gaussian": GaussianModel,
    "lorentzian": LorentzianModel,
    "voigt": VoigtModel,
}


def _make_mixture(n_peaks: int, model_type: str = "gaussian") -> Model:
    """Build a composite model of n_peaks components."""
    model_cls = MODEL_MAP.get(model_type, GaussianModel)
    if model_cls is None:
        raise ValueError(f"Unknown model: {model_type}. Use: {list(MODEL_MAP)}")

    comps = []
    for i in range(n_peaks):
        prefix = f"p{i}_"
        comp = model_cls(prefix=prefix)
        comps.append(comp)

    model = comps[0]
    for comp in comps[1:]:
        model += comp
    return model


def _guess_peak_params(
    x: np.ndarray, y: np.ndarray, n_peaks: int, model_type: str = "gaussian"
) -> dict:
    """Auto-guess initial parameters for n_peaks."""
    from scipy.signal import find_peaks

    y_abs = np.abs(y)
    peaks, props = find_peaks(y_abs, distance=len(x) // (n_peaks + 1))
    if len(peaks) < n_peaks:
        # Force equal spacing
        step = len(x) // (n_peaks + 1)
        peaks = np.array([(i + 1) * step for i in range(n_peaks)], dtype=int)

    # Use the n largest peaks
    peak_heights = y_abs[peaks]
    top_idx = np.argsort(peak_heights)[-n_peaks:][::-1]
    peaks = peaks[top_idx]
    peaks = np.sort(peaks)

    result = {"centers": [], "amplitudes": [], "sigmas": []}
    for i, p_idx in enumerate(peaks):
        result["centers"].append(x[p_idx])
        result["amplitudes"].append(y[p_idx])
        result["sigmas"].append((x[-1] - x[0]) / (n_peaks * 4))

    return result


def fit_peaks(
    x: np.ndarray,
    y: np.ndarray,
    n_peaks: int,
    model_type: str = "gaussian",
    guess_params: dict | None = None,
) -> dict:
    """Fit n_peaks model to data.

    Parameters
    ----------
    x, y : arrays, input data (e.g., Temperature and DTG/DSC).
    n_peaks : number of component peaks.
    model_type : "gaussian", "lorentzian", or "voigt".
    guess_params : optional pre-configured lmfit Parameters.

    Returns
    -------
    dict with:
      - result: lmfit ModelResult
      - components: list of dicts per peak (center, amplitude, sigma, area, fwhm)
      - x, y: input data
      - best_fit: fitted curve
      - r2: pseudo-R²
    """
    model = _make_mixture(n_peaks, model_type)

    if guess_params is None:
        params = model.make_params()
        guess = _guess_peak_params(x, y, n_peaks, model_type)
        for i in range(n_peaks):
            prefix = f"p{i}_"
            # Set center with bounds
            params[f"{prefix}center"].set(
                value=guess["centers"][i], min=x[0], max=x[-1])
            # Set amplitude
            params[f"{prefix}amplitude"].set(value=guess["amplitudes"][i])
            # Set sigma with bounds
            sigma_min = 0.5 if model_type == "lorentzian" else 0.1
            params[f"{prefix}sigma"].set(
                value=guess["sigmas"][i], min=sigma_min, max=x[-1] - x[0])
    else:
        params = guess_params

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = model.fit(y, params, x=x, nan_policy="omit")

    # Extract components
    components = []
    for i in range(n_peaks):
        prefix = f"p{i}_"
        center = result.params.get(f"{prefix}center", None)
        amp = result.params.get(f"{prefix}amplitude", None)
        sigma = result.params.get(f"{prefix}sigma", None)
        height = result.params.get(f"{prefix}height", None)

        c_val = center.value if center else np.nan
        s_val = sigma.value if sigma else np.nan
        amp_val = amp.value if amp else np.nan

        # FWHM and area
        if model_type == "gaussian":
            fwhm = 2.35482 * s_val if not np.isnan(s_val) else np.nan
            area = amp_val * s_val * np.sqrt(2 * np.pi) if not np.isnan(s_val) and not np.isnan(amp_val) else np.nan
        elif model_type == "lorentzian":
            fwhm = 2 * s_val if not np.isnan(s_val) else np.nan
            area = np.pi * amp_val * s_val if not np.isnan(s_val) and not np.isnan(amp_val) else np.nan
        elif model_type == "voigt":
            fwhm = 3.6013 * s_val if not np.isnan(s_val) else np.nan  # approx
            area = amp_val * s_val * np.sqrt(2 * np.pi) if not np.isnan(s_val) and not np.isnan(amp_val) else np.nan
        else:
            fwhm = np.nan
            area = np.nan

        components.append({
            "component": i,
            "center": float(c_val),
            "amplitude": float(amp_val),
            "sigma": float(s_val),
            "fwhm": float(fwhm),
            "area": float(area),
        })

    # Pseudo R²
    ss_res = np.nansum((y - result.best_fit) ** 2)
    ss_tot = np.nansum((y - np.nanmean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return {
        "result": result,
        "components": components,
        "x": x,
        "y": y,
        "best_fit": result.best_fit,
        "r2": float(r2),
    }


def deconvolve(
    df: pd.DataFrame,
    T_range: tuple[float, float],
    n_peaks: int,
    model_type: str = "gaussian",
    signal: str = "DTG_percent_per_C",
) -> dict:
    """Deconvolve a thermal curve (DTG or DSC) within a temperature range.

    Parameters
    ----------
    df : DataFrame with Temp_C and signal columns.
    T_range : (T_low, T_high) for fitting.
    n_peaks : number of component peaks.
    model_type : "gaussian", "lorentzian", or "voigt".
    signal : column name to fit.

    Returns
    -------
    dict from fit_peaks.
    """
    mask = (df["Temp_C"] >= T_range[0]) & (df["Temp_C"] <= T_range[1])
    x = df.loc[mask, "Temp_C"].values
    y = df.loc[mask, signal].values

    if len(x) < 10:
        return {"error": f"Too few points in T_range: {len(x)}"}

    return fit_peaks(x, y, n_peaks, model_type)


def deconvolve_dtg(
    df: pd.DataFrame,
    T_range: tuple[float, float] = (190, 640),
    n_components: int = 4,
    model_type: str = "gaussian",
) -> dict:
    """Deconvolve DTG curve into components."""
    return deconvolve(df, T_range, n_components, model_type, "DTG_percent_per_C")


def deconvolve_dsc(
    df: pd.DataFrame,
    T_range: tuple[float, float] = (190, 640),
    n_components: int = 4,
    model_type: str = "gaussian",
) -> dict:
    """Deconvolve DSC curve into components."""
    return deconvolve(df, T_range, n_components, model_type, "DSC_uW_per_mg")


def component_summary(fit_result: dict) -> pd.DataFrame:
    """Convert fit components to a summary DataFrame."""
    if "components" not in fit_result:
        return pd.DataFrame()
    df = pd.DataFrame(fit_result["components"])
    df["model"] = fit_result.get("r2", np.nan)
    if len(df) > 0:
        total_area = df["area"].sum()
        df["area_fraction"] = df["area"] / total_area if total_area > 0 else np.nan
    return df


# --- Constrained pool-aware deconvolution ---

def deconvolve_pools(
    df: pd.DataFrame,
    pools_scheme: str | list[dict] = "filimonenko2025",
    T_range: tuple[float, float] = (190, 640),
    model_type: str = "gaussian",
    signal: str = "DTG_percent_per_C",
) -> dict:
    """Fit one Gaussian peak independently in each thermal pool range.

    Each pool is fitted separately with a single Gaussian constrained
    to stay within the pool's temperature boundaries. This ensures:
    - Peak 0 always represents labile, Peak 1 = stable, etc.
    - Per-pool R^2 indicates how well a single Gaussian describes that pool
    - Peak center, width, and area are directly comparable across samples

    Parameters
    ----------
    df : DataFrame with Temp_C and signal columns.
    pools_scheme : pool definition name or list of pool dicts.
    T_range : overall temperature range.
    model_type : "gaussian" (symmetric) or "lorentzian" (wider tails).
    signal : column name to fit.

    Returns
    -------
    dict with peaks list, per-pool fit curves, and overall composite.
    """
    from .pools import get_pool_scheme
    from lmfit.models import GaussianModel, LorentzianModel
    from lmfit import Parameters

    pools = get_pool_scheme(pools_scheme)
    model_cls = GaussianModel if model_type == "gaussian" else LorentzianModel

    mask = (df["Temp_C"] >= T_range[0]) & (df["Temp_C"] <= T_range[1])
    x_full = df.loc[mask, "Temp_C"].values.astype(float)
    y_full = df.loc[mask, signal].values.astype(float)

    if len(x_full) < 20:
        return {"error": "Too few data points"}

    peaks = []
    components_y = []
    composite_y = np.zeros_like(y_full)
    total_area = 0.0

    for i, p in enumerate(pools):
        T_low, T_high = p["T_low"], p["T_high"]

        # Extract data within this pool's range
        pool_mask = (x_full >= T_low) & (x_full <= T_high)
        x_pool = x_full[pool_mask]
        y_pool = y_full[pool_mask]

        if len(x_pool) < 10:
            peaks.append({
                "pool": p["pool"], "T_low_C": T_low, "T_high_C": T_high,
                "T_center_C": np.nan, "amplitude": np.nan, "sigma_C": np.nan,
                "fwhm_C": np.nan, "area": np.nan, "r2": np.nan,
                "area_fraction": np.nan, "fit_quality": "too_few_points",
            })
            components_y.append(np.zeros_like(x_full))
            continue

        T_mid = (T_low + T_high) / 2
        T_width = T_high - T_low

        # Amplitude: estimate from max absolute DTG in pool range
        amp_guess = float(np.max(np.abs(y_pool))) if len(y_pool) > 0 else 0.01
        # Check for near-zero signal (flat DTG in this pool)
        signal_range = float(np.max(y_pool) - np.min(y_pool)) if len(y_pool) > 1 else 0

        # Independent single-Gaussian fit for this pool
        model = model_cls(prefix="p_")
        params = Parameters()
        params.add("p_center", value=T_mid, min=T_low, max=T_high)
        params.add("p_sigma", value=T_width / 6, min=T_width / 20, max=T_width / 2)
        params.add("p_amplitude", value=amp_guess, min=0)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = model.fit(y_pool, params, x=x_pool, nan_policy="omit")

        # Extract parameters
        center = result.params.get("p_center")
        sigma = result.params.get("p_sigma")
        amp = result.params.get("p_amplitude")

        c_val = center.value if center else np.nan
        s_val = sigma.value if sigma else np.nan
        a_val = amp.value if amp else np.nan

        if model_type == "gaussian":
            fwhm = 2.35482 * s_val
            area = a_val * s_val * np.sqrt(2 * np.pi)
        else:
            fwhm = 2 * s_val
            area = np.pi * a_val * s_val

        # Per-pool R^2
        ss_res = np.nansum((y_pool - result.best_fit) ** 2)
        ss_tot = np.nansum((y_pool - np.nanmean(y_pool)) ** 2)
        pool_r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

        # Assess fit quality
        if pool_r2 < 0:
            fit_quality = "poor_negative_r2"
        elif pool_r2 < 0.5:
            fit_quality = "marginal"
        elif pool_r2 < 0.8:
            fit_quality = "acceptable"
        else:
            fit_quality = "good"
        if signal_range < 1e-5:
            fit_quality = "flat_signal"

        peaks.append({
            "pool": p["pool"],
            "T_low_C": T_low,
            "T_high_C": T_high,
            "T_center_C": float(c_val),
            "amplitude": float(a_val),
            "sigma_C": float(s_val),
            "fwhm_C": float(fwhm),
            "area": float(area),
            "r2": float(pool_r2) if not np.isnan(pool_r2) else None,
            "fit_quality": fit_quality,
        })

        if not np.isnan(area):
            total_area += area

        # Evaluate this peak on the full x-range for composite
        comp_full = result.eval(x=x_full)
        components_y.append(comp_full)
        composite_y += comp_full

    # Area fractions
    for pe in peaks:
        pe["area_fraction"] = pe["area"] / total_area if total_area > 0 else np.nan

    # Overall R^2 of composite
    ss_res = np.nansum((y_full - composite_y) ** 2)
    ss_tot = np.nansum((y_full - np.nanmean(y_full)) ** 2)
    overall_r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return {
        "peaks": peaks,
        "best_fit": composite_y,
        "r2": float(overall_r2),
        "x": x_full,
        "y": y_full,
        "components_y": components_y,
    }


def pool_deconv_summary(fit_result: dict) -> pd.DataFrame:
    """Convert constrained pool deconvolution result to DataFrame."""
    if "peaks" not in fit_result:
        return pd.DataFrame()
    return pd.DataFrame(fit_result["peaks"])


def batch_deconvolve_pools(
    samples: dict[str, pd.DataFrame],
    pools_scheme: str = "filimonenko2025",
    T_range: tuple[float, float] = (190, 640),
    signal: str = "DTG_percent_per_C",
) -> pd.DataFrame:
    """Run constrained pool deconvolution on all samples.

    Returns long-format DataFrame with one row per pool per sample.
    """
    frames = []
    for name, df in samples.items():
        try:
            result = deconvolve_pools(df, pools_scheme, T_range, signal=signal)
            if "peaks" in result:
                p_df = pd.DataFrame(result["peaks"])
                p_df["sample"] = name
                p_df["overall_r2"] = result.get("r2", np.nan)
                frames.append(p_df)
        except Exception as e:
            warnings.warn(f"Constrained deconv failed for {name}: {e}")

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)
