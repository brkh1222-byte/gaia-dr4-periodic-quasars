"""Gaia helpers: DR3 mean magnitudes, mock DR4 cadence, and frozen DR4 cleaning rules.

BLINDNESS: before the pre-registration commit nothing here touches Gaia epoch photometry.
`dr3_mean_g` queries only gaiadr3.gaia_source (source_id, position, phot_g_mean_mag) for the
G < 20.7 testability cut. It never reads DR3 epoch photometry or any variability table.
`clean_dr4_epochs` is the frozen post-release rule set; it is tested on mock data only.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402

GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap"


# ---------------------------------------------------------------------------
# DR3 mean G (allowed: a mean over 2014.6-2017.4, no time series)
# ---------------------------------------------------------------------------
def dr3_mean_g(cand: pd.DataFrame, radius_arcsec=1.0) -> pd.DataFrame:
    import pyvo
    svc = pyvo.dal.TAPService(GAIA_TAP)
    out = []
    for _, r in cand.iterrows():
        q = ("SELECT source_id, ra, dec, phot_g_mean_mag, "
             f"DISTANCE(POINT({r.ra_use:.7f},{r.dec_use:+.7f}), POINT(ra,dec))*3600 AS sep "
             "FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec), "
             f"CIRCLE({r.ra_use:.7f},{r.dec_use:+.7f},{radius_arcsec/3600:.7f})) ORDER BY sep")
        for k in range(4):
            try:
                t = svc.run_sync(q).to_table().to_pandas()
                break
            except Exception:
                time.sleep(5 * (k + 1))
        else:
            t = None
        rec = {"cand_id": r.cand_id}
        if t is None:
            rec["dr3_status"] = "query_failed"
        elif len(t) == 0:
            rec["dr3_status"] = "no_source_within_1arcsec"
        else:
            rec.update(dr3_source_id=int(t.source_id[0]), dr3_sep_arcsec=float(t.sep[0]),
                       dr3_g_mean=float(t.phot_g_mean_mag[0]) if np.isfinite(t.phot_g_mean_mag[0]) else np.nan,
                       dr3_n_within_1arcsec=len(t), dr3_status="ok")
        out.append(rec)
        time.sleep(0.2)
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# Mock DR4 cadence (planning only; replaced by real DR4 timestamps after release)
# ---------------------------------------------------------------------------
def per_transit_g_error(g):
    """Approximate per-FoV-transit G uncertainty [mag] (planning assumption, not Gaia data).

    0.03 mag at G = 20, scaling as 10^(0.25 (G - 20)), floor 4 mmag.
    """
    return np.maximum(0.004, 0.03 * 10 ** (0.25 * (np.asarray(g, float) - 20.0)))


def mock_cadence(g_mag, ecl_lat_deg, rng, t0=C.GAIA_DR4_START_MJD, t1=C.GAIA_DR4_END_MJD):
    """Return (visit_times, visit_errors, n_transits) for one source.

    Visit groups recur every 30-90 d (scanning-law precession); each group contains 2-6 FoV
    transits within one day (pairs 106.5 min apart, sometimes on consecutive 6-h spins).
    The total number of transits over 66 months is ~70 at the ecliptic and up to ~120 near
    |beta| = 45 deg, as in the Gaia scanning law. After 1-day binning each group is one visit.
    """
    beta = abs(ecl_lat_deg)
    n_target = 70 + 50 * np.exp(-((beta - 45.0) / 15.0) ** 2)
    t, times, ntr = t0 + rng.uniform(0, 60), [], []
    while t < t1:
        k = int(rng.integers(2, 7))
        times.append(t); ntr.append(k)
        t += rng.uniform(30, 90)
    ntr = np.array(ntr, float)
    ntr = np.maximum(1, np.round(ntr * n_target / ntr.sum()))
    err = per_transit_g_error(g_mag) / np.sqrt(ntr)
    return np.array(times), err, ntr.astype(int)


# ---------------------------------------------------------------------------
# Frozen DR4 cleaning rules (applied after release)
# ---------------------------------------------------------------------------
def clean_dr4_epochs(df: pd.DataFrame, win=(None, None)):
    """Per-transit G photometry -> visit-binned (t_mjd, mag, err).

    Expected columns (DR3 naming; mapped on release day if DR4 differs, see prereg):
      g_transit_time [BJD(TCB) - 2455197.5], g_transit_flux, g_transit_flux_error,
      rejected_by_photometry, rejected_by_variability (booleans, optional).
    Rules: drop rejected transits, non-positive or non-finite flux; mag = -2.5 log10(flux)
    + ZP (ZP cancels in the free offset); err = 1.0857 sigma_F / F; convert time to MJD;
    keep transits inside the test window; bin transits closer than VISIT_GAP_D;
    drop visits > GAIA_OUTLIER_SIGMA robust sigma from the median.
    """
    d = df.copy()
    for flag in ("rejected_by_photometry", "rejected_by_variability"):
        if flag in d:
            d = d[~d[flag].astype(bool)]
    f, fe = d["g_transit_flux"].to_numpy(float), d["g_transit_flux_error"].to_numpy(float)
    ok = np.isfinite(f) & np.isfinite(fe) & (f > 0) & (fe > 0)
    t = d["g_transit_time"].to_numpy(float)[ok] + 2455197.5 - 2400000.5
    mag = -2.5 * np.log10(f[ok]) + 25.6874   # DR3 G zero point; any constant cancels
    err = 1.0857 * fe[ok] / f[ok]
    lo, hi = win
    m = np.ones_like(t, bool)
    if lo is not None:
        m &= t >= lo
    if hi is not None:
        m &= t <= hi
    t, mag, err = t[m], mag[m], err[m]
    if len(t) == 0:
        return t, mag, err
    from lightcurves import bin_nights, clip_outliers
    T, Y, E, _ = bin_nights(t, mag, err, gap=C.VISIT_GAP_D)
    return clip_outliers(T, Y, E, nsig=C.GAIA_OUTLIER_SIGMA)
