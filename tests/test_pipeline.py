"""Tests for light-curve handling, the ephemeris summary and the frozen Gaia cleaning rules.
Mock data only: no Gaia data are read."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import models as M  # noqa: E402
import lightcurves as LCm  # noqa: E402
import gaia as G  # noqa: E402
import frozen as F  # noqa: E402


def test_bin_nights_weighted_mean():
    t = np.array([1.0, 1.1, 1.2, 5.0])
    y = np.array([10.0, 10.2, 10.1, 11.0])
    e = np.array([0.1, 0.1, 0.1, 0.2])
    T, Y, E, N = LCm.bin_nights(t, y, e)
    assert len(T) == 2 and N.tolist() == [3, 1]
    assert np.isclose(Y[0], 10.1) and np.isclose(Y[1], 11.0)
    assert E[0] >= 0.1 / np.sqrt(3) - 1e-12


def test_clip_outliers_removes_spike():
    rng = np.random.default_rng(0)
    t = np.arange(100.0); y = 18 + 0.05 * rng.standard_normal(100); y[50] = 25
    e = np.full(100, 0.05)
    t2, y2, _ = LCm.clip_outliers(t, y, e)
    assert 50 not in t2 and len(t2) >= 97


def test_ephemeris_recovers_brightest_epoch():
    """Noise-free strong sinusoid: the brightest epoch must match the truth to ~1 d."""
    P, t_bright_true = 1000.0, 55123.0
    t = np.linspace(53000, 56500, 400)
    y = 18 + 0.3 * np.cos(2 * np.pi * (t - t_bright_true) / P + np.pi)   # min mag at t_bright
    lc = M.LightCurve(t, y, np.full_like(t, 0.01))
    hS = M.Hyper(sigma=0.01, tau=50.0, jitter=0.0, sigma_s=0.3, period=P)
    amp, amp_sd, tb, tb_sd = F.ephemeris(lc, hS, 57000.0)
    k = np.round((tb - t_bright_true) / P)
    assert abs(tb - t_bright_true - k * P) < 2.0 and tb >= 57000.0 and tb < 57000.0 + P
    assert abs(amp - 0.3) < 0.01


def test_mock_cadence_reasonable():
    rng = np.random.default_rng(1)
    t, e, ntr = G.mock_cadence(19.5, 30.0, rng)
    assert 20 <= len(t) <= 50 and 60 <= ntr.sum() <= 140
    assert np.all(np.diff(t) > 1.0) and e.max() < 0.05


def test_clean_dr4_epochs_on_mock_table():
    rng = np.random.default_rng(2)
    n = 90
    tt = np.sort(rng.uniform(1700, 3700, n))          # BJD - 2455197.5
    flux = 10 ** (-0.4 * (19.0 - 25.6874)) * (1 + 0.02 * rng.standard_normal(n))
    df = pd.DataFrame(dict(g_transit_time=tt, g_transit_flux=flux,
                           g_transit_flux_error=0.02 * flux,
                           rejected_by_photometry=False, rejected_by_variability=False))
    df.loc[3, "rejected_by_photometry"] = True
    df.loc[5, "g_transit_flux"] = -1.0
    T, Y, E = G.clean_dr4_epochs(df)
    assert np.all(np.diff(T) > 1.0)
    assert abs(np.median(Y) - 19.0) < 0.02
    assert T.min() > 56863 - 1 and T.max() < 58870
