"""Fit and freeze the per-candidate model parameters (task T2.4 input).

For every candidate in the scored groups (forecast: Graham+2015, Charisi+2016; gap:
Foustoul+2025) that has usable discovery data and Gaia DR3 mean G < 20.7:

  * model D (DRW) and model S (sinusoid at the published period + DRW) hyperparameters are
    fitted to the discovery data only (MAP, src/models.py) and frozen;
  * the sinusoid's phase and amplitude posterior (given the discovery data only) is
    summarised as an ephemeris: semi-amplitude A and the MJD of maximum brightness
    (minimum magnitude) nearest to the start of the test window;
  * the frozen predictive mean and sd of both models are tabulated every 10 d across the
    test window (nuisances at s = 1, c = 0; in the discovery band);
  * the expected power is estimated on a mock Gaia cadence: the fraction of model-S
    simulated data sets whose conditional p-value is < 0.05, and P(lnB > 3).

Outputs (written to results/ by default; copied to preregistration/ for the commit):
  frozen_params.csv, frozen_predictions.csv, power_mock.csv, discovery_summary.csv

Usage: .venv/bin/python src/frozen.py [--outdir results/frozen] [--n-sim 1000]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, GeocentricTrueEcliptic
import astropy.units as u

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import models as M  # noqa: E402
import score as S  # noqa: E402
import lightcurves as LCm  # noqa: E402
import gaia as G  # noqa: E402

SCORED_SAMPLES = ("Graham2015", "Charisi2016", "Foustoul2025")
DEFAULT_OUT = C.RESULTS / "frozen"


def cid_seed(cid: str) -> int:
    return int(hashlib.sha256(cid.encode()).hexdigest()[:8], 16)


def ecl_lat(ra, dec):
    c = SkyCoord(ra * u.deg, dec * u.deg).transform_to(GeocentricTrueEcliptic())
    return float(c.lat.deg)


def mock_gaia(row, win):
    rng = np.random.default_rng(cid_seed(row["cand_id"]))
    t, e, ntr = G.mock_cadence(row["dr3_g_mean"], ecl_lat(row["ra_deg"], row["dec_deg"]), rng)
    m = (t >= win[0]) & (t <= win[1])
    return t[m], e[m]


def ephemeris(lc, hypS, t_ref):
    """Semi-amplitude and epoch of maximum brightness from the posterior mean of (a, b)."""
    cd = M.condition(lc, hypS, np.array([t_ref]), hypS.period)
    ns = lc.n_surveys
    a, b = cd.m[ns], cd.m[ns + 1]
    S_ab = cd.S[ns:ns + 2, ns:ns + 2]
    amp = float(np.hypot(a, b))
    # x_sin(t) = a cos(wt) + b sin(wt) = amp cos(wt - phi), phi = atan2(b, a).
    # Magnitude is minimal (brightest) where cos(wt - phi) = -1.
    w = 2 * np.pi / hypS.period
    phi = np.arctan2(b, a)
    t_bright = (phi + np.pi) / w
    k = np.ceil((t_ref - t_bright) / hypS.period)
    t_bright += k * hypS.period
    # amplitude uncertainty (delta method) and phase uncertainty in days
    g = np.array([a, b]) / max(amp, 1e-12)
    amp_sd = float(np.sqrt(g @ S_ab @ g))
    gp = np.array([-b, a]) / max(amp**2, 1e-24)
    phase_sd_d = float(np.sqrt(gp @ S_ab @ gp) / w)
    return amp, amp_sd, float(t_bright), phase_sd_d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(DEFAULT_OUT))
    ap.add_argument("--n-sim", type=int, default=1000)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    cand = pd.read_csv(C.CANDIDATES_CSV)
    dr3 = pd.read_csv(C.RAW / "dr3_mean_g.csv")
    cand = cand.merge(dr3, on="cand_id", how="left")
    if a.only:
        cand = cand[cand.cand_id.isin(a.only.split(","))]
    params, preds, power, summ = [], [], [], []
    for _, row in cand.iterrows():
        rec = {"cand_id": row.cand_id, "sample": row["sample"], "group": row["group"],
               "name": row["name"], "dr3_g_mean": row.get("dr3_g_mean", np.nan)}
        if row["sample"] not in SCORED_SAMPLES:
            summ.append({**rec, "status": "not scored (exploratory sample)"}); continue
        if row["sample"] == "Charisi2016" and not bool(row["in_final_sample"]):
            summ.append({**rec, "status": "not scored (not in Charisi final 33)"}); continue
        if not np.isfinite(row["period_d"]):
            summ.append({**rec, "status": "no published period"}); continue
        g = row.get("dr3_g_mean", np.nan)
        if not (np.isfinite(g) and g < C.GAIA_G_LIMIT):
            summ.append({**rec, "status": f"G cut (DR3 mean G = {g})"}); continue
        lc, win, info = LCm.discovery_data(row)
        if lc is None:
            summ.append({**rec, "status": info.get("reason", "no discovery data")}); continue
        P = float(row["period_d"])
        sp = row["period_err_d"] if np.isfinite(row["period_err_d"]) else C.SIGMA_P_FRAC_DEFAULT * P
        rng = np.random.default_rng(cid_seed(row.cand_id) + 1)
        hD = M.fit_hyper(lc, "D", rng=rng)
        hS = M.fit_hyper(lc, "S", period=P, sigma_p=sp, rng=rng)
        lzS = M.log_evidence(lc, hS.sigma, hS.tau, hS.jitter, P, hS.sigma_s)
        lzD = M.log_evidence(lc, hD.sigma, hD.tau, hD.jitter)
        amp, amp_sd, t_bright, ph_sd = ephemeris(lc, hS, win[0])
        params.append(dict(
            cand_id=row.cand_id, sample=row["sample"], group=row["group"], name=row["name"],
            ra_deg=row.ra_deg, dec_deg=row.dec_deg, period_d=P, period_err_sensitivity_d=sp,
            disc_survey=info["survey"], n_disc_nights=len(lc.t), n_surveys=lc.n_surveys,
            disc_t_first=lc.t.min(), disc_t_last=lc.t.max(),
            test_win_lo=win[0], test_win_hi=win[1],
            D_sigma=hD.sigma, D_tau=hD.tau, D_jitter=hD.jitter,
            S_sigma=hS.sigma, S_tau=hS.tau, S_jitter=hS.jitter, S_sigma_s=hS.sigma_s,
            sin_semi_amp=amp, sin_semi_amp_sd=amp_sd, sin_t_brightest_mjd=t_bright,
            sin_t_brightest_sd_d=ph_sd, lnZ_S_minus_lnZ_D_disc=lzS - lzD,
            D_fit_converged=hD.extra["converged"], S_fit_converged=hS.extra["converged"],
            dr3_g_mean=g))
        # frozen predictive curves
        tg = np.arange(win[0], win[1] + 1e-6, 10.0)
        mS, sS = M.predictive_moments(lc, hS, tg, P)
        mD, sD = M.predictive_moments(lc, hD, tg, None)
        preds.append(pd.DataFrame(dict(cand_id=row.cand_id, mjd=tg, S_mean=mS, S_sd=sS,
                                       D_mean=mD, D_sd=sD)))
        # expected power on the mock cadence
        t_g, e_g = mock_gaia(row, win)
        if len(t_g) >= 5:
            sc = S.CandidateScorer(lc, hS, hD, t_g, e_g)
            nul = sc.null_lnB(a.n_sim, rng)
            inj = sc.inj_lnB(a.n_sim, rng)
            crit = np.quantile(nul, 0.95)
            power.append(dict(cand_id=row.cand_id, sample=row["sample"], n_mock_visits=len(t_g),
                              null_q95=crit, power_p05=float(np.mean(inj > crit)),
                              power_lnB3=float(np.mean(inj > C.LNB_THRESHOLD)),
                              null_frac_lnB3=float(np.mean(nul > C.LNB_THRESHOLD)),
                              inj_median_lnB=float(np.median(inj)),
                              null_median_lnB=float(np.median(nul))))
        summ.append({**rec, "status": "scored", "n_disc_nights": len(lc.t),
                     "test_window_d": win[1] - win[0]})
        print(f"{row.cand_id} {row['name'][:28]:28s} G={g:.2f} n={len(lc.t)} "
              f"sig_s={hS.sigma_s:.3f} A={amp:.3f}±{amp_sd:.3f} dlnZ={lzS-lzD:+.2f} "
              f"pow={power[-1]['power_p05'] if power and power[-1]['cand_id']==row.cand_id else np.nan:.2f}",
              flush=True)
    pd.DataFrame(params).to_csv(out / "frozen_params.csv", index=False, float_format="%.6g")
    if preds:
        pd.concat(preds).to_csv(out / "frozen_predictions.csv", index=False, float_format="%.5f")
    pd.DataFrame(power).to_csv(out / "power_mock.csv", index=False, float_format="%.4g")
    pd.DataFrame(summ).to_csv(out / "discovery_summary.csv", index=False)
    print(pd.DataFrame(summ).groupby(["sample", "status"]).size())


# ---------------------------------------------------------------------------
# Template access for simulate_null.py
# ---------------------------------------------------------------------------
def _frozen_table(outdir=DEFAULT_OUT):
    return pd.read_csv(Path(outdir) / "frozen_params.csv")


def template_ids(outdir=DEFAULT_OUT, groups=("forecast",)):
    fp = _frozen_table(outdir)
    return fp[fp.group.isin(groups)].cand_id.tolist()


def load_template(cid, outdir=DEFAULT_OUT):
    fp = _frozen_table(outdir).set_index("cand_id")
    p = fp.loc[cid]
    cand = pd.read_csv(C.CANDIDATES_CSV).merge(pd.read_csv(C.RAW / "dr3_mean_g.csv"), on="cand_id")
    row = cand.set_index("cand_id").loc[cid].copy()
    row["cand_id"] = cid
    lc, win, info = LCm.discovery_data(row)
    hD = M.Hyper(sigma=p.D_sigma, tau=p.D_tau, jitter=p.D_jitter)
    hS = M.Hyper(sigma=p.S_sigma, tau=p.S_tau, jitter=p.S_jitter, sigma_s=p.S_sigma_s,
                 period=p.period_d, sigma_p=p.period_err_sensitivity_d)
    t_g, e_g = mock_gaia(row, win)
    return row, lc, hS, hD, t_g, e_g


if __name__ == "__main__":
    main()
