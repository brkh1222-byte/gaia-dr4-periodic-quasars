"""Pre-registered DR4 analysis (tasks T3.3; frozen with the pre-registration).

Input : data/gaia_epochs/<cand_id>.csv  per-transit G photometry (DR3-style columns, see
        gaia.clean_dr4_epochs), produced by src/gaia_epochs.py after release
        results/frozen/frozen_params.csv (frozen hyperparameters; identical to the copy in
        preregistration/)
Output: results/<tag>/scores.csv, results/<tag>/population.json

Steps (all rules from preregistration/prereg.md):
  1. clean and visit-bin the Gaia epochs in each candidate's test window;
  2. testability from timestamps and errors only: power = fraction of 1000 model-S draws
     whose conditional null p < 0.05 (the Gaia magnitudes are not used in this step);
  3. ln B for the observed Gaia data; conditional null p from 1000 model-D draws;
  4. H1 binomial test, H2 (ln B > 3 and BH q < 0.05), H3 posterior on f, H4 phase scrambling
     for survivors; controls and sensitivity analyses on request.

`--mock null|signal` runs the dress rehearsal: it writes mock epoch files drawn from model D
or model S on the mock cadence (no Gaia data involved) and analyses them.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import frozen as F  # noqa: E402
import gaia as G  # noqa: E402
import lightcurves as LCm  # noqa: E402
import models as M  # noqa: E402
import score as S  # noqa: E402

TESTABLE_POWER_P05 = 0.5   # frozen testability rule (prereg.md §2; notes/decisions.md)


def load_frozen(fp_row):
    hD = M.Hyper(sigma=fp_row.D_sigma, tau=fp_row.D_tau, jitter=fp_row.D_jitter)
    hS = M.Hyper(sigma=fp_row.S_sigma, tau=fp_row.S_tau, jitter=fp_row.S_jitter,
                 sigma_s=fp_row.S_sigma_s, period=fp_row.period_d,
                 sigma_p=fp_row.period_err_sensitivity_d)
    return hS, hD


def write_mock_epochs(fp, outdir, kind, seed=12345):
    """Mock per-transit files in DR3 format from model D or S predictive (dress rehearsal)."""
    outdir.mkdir(parents=True, exist_ok=True)
    cand = pd.read_csv(C.CANDIDATES_CSV).merge(pd.read_csv(C.RAW / "dr3_mean_g.csv"), on="cand_id")
    cand = cand.set_index("cand_id")
    for _, p in fp.iterrows():
        row = cand.loc[p.cand_id].copy(); row["cand_id"] = p.cand_id
        lc, win, _ = LCm.discovery_data(row)
        hS, hD = load_frozen(p)
        t_g, e_g = F.mock_gaia(row, win)
        if len(t_g) < 5:
            continue
        sc = S.CandidateScorer(lc, hS, hD, t_g, e_g)
        rng = np.random.default_rng(seed + F.cid_seed(p.cand_id) % 100000)
        y = (sc.D if kind == "null" else sc.S).sample(1, rng)[0]
        # one transit per visit with the visit error (binning then returns the same values)
        flux = 10 ** (-0.4 * (y - 25.6874))
        df = pd.DataFrame(dict(g_transit_time=t_g - (2455197.5 - 2400000.5),
                               g_transit_flux=flux, g_transit_flux_error=flux * e_g / 1.0857,
                               rejected_by_photometry=False, rejected_by_variability=False))
        df.to_csv(outdir / f"{p.cand_id}.csv", index=False)


SENSITIVITY = {
    "primary": {},
    "dr4only": {"dr4_only": True},            # epochs after the end of DR3 (2017-05-28)
    "unbinned": {"visit_gap": 0.0},           # individual transits, no visit binning
    "pmarg": {"marginalise_period": True},    # period marginalised over +-3 sigma_P
    "widescale": {"scale_sigma_ln": 0.7},     # wider G/survey amplitude-ratio prior
    "shuffle": {"shuffle_periods": True},     # negative control: periods permuted
}


def analyse(fp, epoch_dir, outdir, n_sim=1000, groups=("forecast",), seed=C.RANDOM_SEED,
            dr4_only=False, h4=True, visit_gap=None, marginalise_period=False,
            scale_sigma_ln=None, shuffle_periods=False):
    outdir.mkdir(parents=True, exist_ok=True)
    cand = pd.read_csv(C.CANDIDATES_CSV).merge(pd.read_csv(C.RAW / "dr3_mean_g.csv"), on="cand_id")
    cand = cand.set_index("cand_id")
    rows, nulls, injs = [], {}, {}
    fp = fp[fp.group.isin(groups)].copy()
    if shuffle_periods:   # negative control: each candidate gets another candidate's period
        perm = np.random.default_rng(seed + 7).permutation(len(fp))
        fp["period_d"] = fp["period_d"].values[perm]
        fp["period_err_sensitivity_d"] = fp["period_err_sensitivity_d"].values[perm]
    for _, p in fp.iterrows():
        rec = dict(cand_id=p.cand_id, sample=p["sample"], group=p.group, period_d=p.period_d)
        f = epoch_dir / f"{p.cand_id}.csv"
        if not f.exists():
            rows.append({**rec, "status": "no DR4 epochs"}); continue
        row = cand.loc[p.cand_id].copy(); row["cand_id"] = p.cand_id
        lc, win, _ = LCm.discovery_data(row)
        lo = max(win[0], C.GAIA_DR3_END_MJD) if dr4_only else win[0]
        t_g, y_g, e_g = G.clean_dr4_epochs(pd.read_csv(f), (lo, win[1]), visit_gap=visit_gap)
        rec["n_visits"] = len(t_g)
        if len(t_g) < 5:
            rows.append({**rec, "status": "fewer than 5 visits"}); continue
        hS, hD = load_frozen(p)
        if shuffle_periods:   # the S-model hyperparameters are refitted at the swapped period
            hS = M.fit_hyper(lc, "S", period=p.period_d, sigma_p=p.period_err_sensitivity_d,
                             rng=np.random.default_rng(F.cid_seed(p.cand_id) + 1))
        sc = S.CandidateScorer(lc, hS, hD, t_g, e_g, marginalise_period=marginalise_period,
                               scale_sigma_ln=scale_sigma_ln)
        rng = np.random.default_rng(seed + F.cid_seed(p.cand_id) % 100000)
        # (2) testability: timestamps and errors only
        nul = sc.null_lnB(n_sim, rng)
        inj = sc.inj_lnB(n_sim, rng)
        crit = np.quantile(nul, 0.95)
        power = float(np.mean(inj > crit))
        # (3) score
        lnb = float(sc.lnB(y_g[None])[0])
        pval = S.null_pvalue(lnb, nul)
        rows.append({**rec, "status": "scored", "power_p05": power,
                     "testable": power >= TESTABLE_POWER_P05, "lnB": lnb, "p": pval,
                     "power_lnB3": float(np.mean(inj > C.LNB_THRESHOLD))})
        nulls[p.cand_id], injs[p.cand_id] = nul, inj
    df = pd.DataFrame(rows)
    pop = {}
    T = df[(df.get("status") == "scored") & (df.get("testable") == True)] if len(df) else df  # noqa: E712
    if len(T):
        r0_used = json.load(open(C.RESULTS / "null_sim_summary.json"))["r0_used"] \
            if (C.RESULTS / "null_sim_summary.json").exists() else C.P_SINGLE
        n, k = len(T), int((T.p < C.P_SINGLE).sum())
        pop["H1"] = dict(n_testable=n, k=k, rate=r0_used, expected=n * r0_used,
                         p_binom=S.binomial_excess_pvalue(k, n, r0_used),
                         reject_H0=S.binomial_excess_pvalue(k, n, r0_used) < C.ALPHA_POP)
        bh = S.benjamini_hochberg(T.p.values, C.FDR_Q)
        surv = T[(T.lnB > C.LNB_THRESHOLD).values & bh]
        df["bh_discovery"] = df.cand_id.isin(T.cand_id[bh])
        h4res = {}
        if h4:
            for cid in surv.cand_id:
                p = fp.set_index("cand_id").loc[cid]
                row = cand.loc[cid].copy(); row["cand_id"] = cid
                lc, win, _ = LCm.discovery_data(row)
                t_g, y_g, e_g = G.clean_dr4_epochs(pd.read_csv(epoch_dir / f"{cid}.csv"), win)
                hS, hD = load_frozen(p)
                rng = np.random.default_rng(seed + 99)
                phis = rng.uniform(0, 2 * np.pi, C.N_PHASE_SCRAMBLES)
                base = S.CandidateScorer(lc, hS, hD, t_g, e_g).lnB(y_g[None])[0]
                scr = [S.CandidateScorer(lc, hS, hD, t_g, e_g, phase_rotation=ph).lnB(y_g[None])[0]
                       for ph in phis]
                per = {fac: float(S.CandidateScorer(lc, hS, hD, t_g, e_g, period=fac * hS.period)
                                  .lnB(y_g[None])[0]) for fac in (0.9, 1.1)}
                h4res[cid] = dict(lnB=float(base), frac_scrambled_below=float(np.mean(np.array(scr) < base)),
                                  passes=bool(np.mean(np.array(scr) < base) >= 0.95), lnB_period=per)
        pop["H2"] = dict(n_lnB_gt3=int((T.lnB > C.LNB_THRESHOLD).sum()), n_bh=int(bh.sum()),
                         survivors=list(surv.cand_id),
                         ephemeris_specific=[c for c in surv.cand_id if h4res.get(c, {}).get("passes")],
                         h4=h4res)
        fpost = S.fraction_posterior(T.lnB.values, [nulls[c] for c in T.cand_id],
                                     [injs[c] for c in T.cand_id])
        pop["H3"] = dict(f_median=fpost["median"], f_upper95=fpost["upper95"])
    df.to_csv(outdir / "scores.csv", index=False)
    json.dump(pop, open(outdir / "population.json", "w"), indent=2, default=float)
    return df, pop


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mock", choices=["null", "signal"], default=None)
    ap.add_argument("--n-sim", type=int, default=1000)
    ap.add_argument("--tag", default="dr4")
    ap.add_argument("--groups", default="forecast")
    ap.add_argument("--dr4-only", action="store_true")
    ap.add_argument("--sensitivity", default="primary", choices=list(SENSITIVITY))
    a = ap.parse_args()
    fp = pd.read_csv(C.RESULTS / "frozen" / "frozen_params.csv")
    if a.mock:
        epoch_dir = C.DATA / f"mock_epochs_{a.mock}"
        write_mock_epochs(fp, epoch_dir, a.mock)
        tag = f"dress_rehearsal_{a.mock}"
    else:
        epoch_dir, tag = C.DATA / "gaia_epochs", a.tag
    opts = dict(SENSITIVITY[a.sensitivity])
    if a.dr4_only:
        opts["dr4_only"] = True
    if a.sensitivity != "primary":
        tag = f"{tag}_{a.sensitivity}"
    df, pop = analyse(fp, epoch_dir, C.RESULTS / tag, a.n_sim, tuple(a.groups.split(",")),
                      h4=(a.sensitivity == "primary"), **opts)
    print(df.groupby(["sample", "status"]).size())
    print(json.dumps(pop, indent=1, default=float)[:3000])


if __name__ == "__main__":
    main()
