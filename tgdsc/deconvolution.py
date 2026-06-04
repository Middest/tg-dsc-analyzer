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
