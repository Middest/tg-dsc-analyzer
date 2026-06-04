"""Unit tests for tg-dsc-analyzer core modules."""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import pytest

from tgdsc.io import read_netzsch_excel, batch_read, merge_samples
from tgdsc.preprocess import smooth_savgol, remove_duplicate_time
from tgdsc.tg import mass_loss, T50, char_temperatures, residual_mass_percent
from tgdsc.dsc import energy_density, detect_peaks, peak_integral
from tgdsc.pools import (
    POOL_SCHEMES, get_pool_scheme, calculate_pools, pool_mass_loss,
)
from tgdsc.kinetics import coats_redfern, pool_kinetics, weighted_ea
from tgdsc.stability import thermostability_index, comprehensive_report
from tgdsc.deconvolution import fit_peaks, deconvolve_dtg, deconvolve_pools, batch_deconvolve_pools


# Fixtures

@pytest.fixture
def sample_data() -> pd.DataFrame:
    """Generate synthetic TG-DSC data for testing."""
    np.random.seed(42)
    n = 500
    T = np.linspace(30, 800, n)
    # Simulate typical TG curve: slow loss then rapid degradation
    TG = 100 - 0.5 * np.exp(-((T - 400) ** 2) / 5000) - 5 * (T / 800) ** 3
    TG += np.random.normal(0, 0.02, n)
    # DTG: derivative of mass loss
    DTG = np.gradient(-TG, T[1] - T[0])
    # DSC: simulate exothermic peaks
    DSC = 2 * np.exp(-((T - 250) ** 2) / 2000) + 5 * np.exp(-((T - 450) ** 2) / 3000) - 5
    DSC += np.random.normal(0, 0.3, n)

    return pd.DataFrame({
        "Time_min": np.linspace(0, 77, n),
        "Temp_C": T,
        "TG_percent": TG,
        "DTG_percent_per_C": DTG,
        "DSC_uW_per_mg": DSC,
    })


@pytest.fixture
def sample_file(sample_data, tmp_path):
    """Create a temporary NETZSCH-format Excel file."""
    # Build the 9-column format
    rows = [["Time", "Temp", "TG", "Time", "Temp", "DTG", "Time", "Temp", "DSC"],
            ["min", "C", "%", "min", "C", "%/Cel", "min", "C", "uW/mg"]]
    for _, row in sample_data.iterrows():
        t = row["Time_min"]
        T = row["Temp_C"]
        tg = row["TG_percent"]
        dtg = row["DTG_percent_per_C"]
        dsc = row["DSC_uW_per_mg"]
        rows.append([t, T, tg, t, T, dtg, t, T, dsc])

    df = pd.DataFrame(rows[1:], columns=rows[0])
    fpath = tmp_path / "test_sample.xlsx"
    df.to_excel(fpath, index=False, header=False)
    return fpath


# IO tests

def test_read_netzsch_excel(sample_file):
    df = read_netzsch_excel(sample_file)
    assert "Temp_C" in df.columns
    assert "TG_percent" in df.columns
    assert len(df) > 100
    assert df["TG_percent"].values[0] > 90


def test_merge_samples(sample_data):
    samples = {"A": sample_data.copy(), "B": sample_data.copy()}
    tmap = {"A": {"Treatment": "CK"}, "B": {"Treatment": "BC"}}
    merged = merge_samples(samples, tmap)
    assert len(merged) == len(sample_data) * 2
    assert "Treatment" in merged.columns


# Preprocessing tests

def test_smooth_savgol(sample_data):
    smoothed = smooth_savgol(sample_data, window=15)
    assert smoothed["TG_percent"].std() < sample_data["TG_percent"].std()


def test_remove_duplicate_time(sample_data):
    dup = pd.concat([sample_data, sample_data.iloc[:10]])
    clean = remove_duplicate_time(dup)
    assert len(clean) == len(sample_data)


# TG tests

def test_mass_loss(sample_data):
    ml = mass_loss(sample_data, 200, 600)
    assert 0 < ml < 50  # Typical SOM mass loss range


def test_T50(sample_data):
    t50 = T50(sample_data, (200, 600))
    assert 200 < t50 < 600


# DSC tests

def test_peak_integral(sample_data):
    result = peak_integral(sample_data, 200, 600)
    assert "energy_J_g_soil" in result


def test_energy_density(sample_data):
    ed = energy_density(sample_data, (200, 600))
    assert "ED_kJ_g_OM" in ed


# Pools tests

def test_get_pool_scheme():
    pools = get_pool_scheme("filimonenko2025")
    assert len(pools) == 4
    assert pools[0]["pool"] == "labile"


def test_calculate_pools(sample_data):
    df = calculate_pools(sample_data, "filimonenko2025", (200, 600))
    assert len(df) == 4
    assert "mass_loss_percent" in df.columns


# Kinetics tests

def test_coats_redfern(sample_data):
    cr = coats_redfern(sample_data, (300, 500))
    assert "Ea_kJ_mol" in cr
    assert "R2" in cr


def test_pool_kinetics(sample_data):
    kin = pool_kinetics(sample_data)
    assert len(kin) == 4


# Stability tests

def test_thermostability_index(sample_data):
    ti = thermostability_index(sample_data)
    assert 0 < ti["TI"] < 1


def test_comprehensive_report(sample_data):
    report = comprehensive_report(sample_data)
    assert "T50_C" in report
    assert "ED_kJ_g_OM" in report


# Deconvolution tests

def test_fit_peaks(sample_data):
    x = sample_data["Temp_C"].values[100:400]
    y = sample_data["DTG_percent_per_C"].values[100:400]
    result = fit_peaks(x, y, n_peaks=3, model_type="gaussian")
    assert "result" in result
    assert "components" in result
    assert len(result["components"]) == 3


def test_deconvolve_pools(sample_data):
    result = deconvolve_pools(sample_data, "filimonenko2025", (200, 600))
    assert "peaks" in result
    assert len(result["peaks"]) == 4
    assert result["peaks"][0]["pool"] == "labile"
    assert "T_center_C" in result["peaks"][0]
    # Each peak center should be within its pool range
    for peak in result["peaks"]:
        if not np.isnan(peak["T_center_C"]):
            assert peak["T_low_C"] <= peak["T_center_C"] <= peak["T_high_C"]


def test_batch_deconvolve_pools(sample_data):
    samples = {"A": sample_data, "B": sample_data.copy()}
    df = batch_deconvolve_pools(samples, "filimonenko2025", (200, 600))
    assert len(df) == 8  # 2 samples x 4 pools
    assert "sample" in df.columns
    assert "pool" in df.columns
    assert "T_center_C" in df.columns
    assert df["sample"].nunique() == 2
