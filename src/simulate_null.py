"""End-to-end DRW null simulation with "fake candidates" (task T2.3).

Purpose: calibrate the whole pipeline (hyperparameter fits, scoring, conditional p-values)
on red-noise light curves that were *selected to look periodic*, as real candidates were.
This yields, before any DR4 data exist,
  * r0  = the false-positive rate P(p < 0.05) per fake candidate (null rate used by the
          H1 binomial test: r0_used = max(0.05, upper 95% bound of r0));
  * the false-survivor rate P(lnB > 3);
  * the distribution of lnB under the selected-DRW null.

For each template (a real candidate's discovery cadence, errors, frozen DRW hyperparameters
and published period P_i) and a mock Gaia DR4 cadence in the candidate's test window:
  1. simulate a DRW light curve at discovery and Gaia epochs (true hyperparameters = the
     candidate's model-D fit); Gaia data get nuisance draws s ~ LogN(0,0.35),
     j ~ HalfN(0.05), c ~ N(0, 0.3^2);
  2. selection: GP periodogram  dlnZ(P) = ln Z_S(P) - ln Z_D  on 41 log-spaced periods in
     [0.8 P_i, 1.25 P_i]; accept if max dlnZ >= dlnZ of the real candidate at P_i
     (the fake looks at least as periodic as the real one, at a similar period).
     If fewer than K fakes are accepted in N_TRIALS_MAX trials, the K with the largest
     dlnZ are used and the template is flagged 'relaxed';
  3. the fake's 'published period' is its periodogram peak P_f; hyperparameters are refitted
     from scratch exactly as for real candidates; lnB and the conditional p-value (N_COND
     draws from the fitted model-D predictive) are computed exactly as for real data.
Optionally (--psd dho) the fakes are generated from a steeper damped-harmonic-oscillator
(CARMA(2,1)) process instead of a DRW, for the pre-registered sensitivity analysis.

Usage:
  .venv/bin/python src/simulate_null.py --k 50 --workers 2 [--templates CH001,GR005] [--psd dho]
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import models as M  # noqa: E402
import score as S  # noqa: E402

N_GRID = 41
GRID_LO, GRID_HI = 0.8, 1.25
N_TRIALS_MAX = 20000
BATCH = 500
N_COND = 200
OFFSET_SIM_SIGMA = 0.3


# ---------------------------------------------------------------------------
# Batched GP periodogram via Woodbury identity
# ---------------------------------------------------------------------------
class Periodogram:
    def __init__(self, t, e, survey, n_surveys, hypS, hypD, periods):
        self.t, self.periods = t, periods
        self.ns = n_surveys
        self.facS = M.Factor(t, e, hypS.sigma, hypS.tau, hypS.jitter)
        self.facD = M.Factor(t, e, hypD.sigma, hypD.tau, hypD.jitter)
        self.HD = M.basis(t, survey, n_surveys)
        self.lamD = M.prior_cov_diag(n_surveys)
        self.HS = [M.basis(t, survey, n_surveys, P) for P in periods]
        self.lamS = M.prior_cov_diag(n_surveys, hypS.sigma_s)
        self._pre()

    @staticmethod
    def _parts(fac, H, lam):
        KiH = fac.solve(H)
        A = np.diag(1.0 / lam) + H.T @ KiH
        Ainv = np.linalg.inv(A)
        sign, ld = np.linalg.slogdet(A)
        const = fac.logdet + ld + np.sum(np.log(lam))
        return H, Ainv, const

    def _pre(self):
        self.pD = self._parts(self.facD, self.HD, self.lamD)
        self.pS = [self._parts(self.facS, H, self.lamS) for H in self.HS]

    def dlnZ(self, Y):
        """Y: (n, B) -> (n_periods, B) array of ln Z_S(P) - ln Z_D."""
        def lnz(fac, parts, KiY=None):
            H, Ainv, const = parts
            KiY = fac.solve(Y) if KiY is None else KiY
            U = H.T @ KiY
            quad = np.sum(Y * KiY, axis=0) - np.sum(U * (Ainv @ U), axis=0)
            return -0.5 * (quad + const), KiY
        zD, _ = lnz(self.facD, self.pD)
        KiYS = self.facS.solve(Y)
        out = np.empty((len(self.periods), Y.shape[1]))
        for k, parts in enumerate(self.pS):
            out[k], _ = lnz(self.facS, parts, KiYS)
        return out - zD[None, :]


# ---------------------------------------------------------------------------
def simulate_dho(t, sigma, tau, rng, n_draws, q=0.4):
    """Steeper-PSD alternative: a damped harmonic oscillator (CARMA(2,1)-like) with the same
    variance sigma^2 and characteristic timescale tau (celerite2 SHOTerm, quality factor q)."""
    import celerite2
    from celerite2 import terms
    w0 = 1.0 / tau
    k = terms.SHOTerm(sigma=sigma, w0=w0, Q=q)
    gp = celerite2.GaussianProcess(k, mean=0.0)
    gp.compute(np.sort(t), diag=np.full(len(t), 1e-10))
    o = np.argsort(t)
    draws = gp.sample(size=n_draws, random_state=int(rng.integers(2**31)))
    out = np.empty((n_draws, len(t)))
    out[:, o] = np.atleast_2d(draws)
    return out


def run_template(row, lc, hypS, hypD, t_g, e_g, k, rng, psd="drw"):
    """Generate k accepted fakes for one template and score them end to end."""
    t_all = np.concatenate([lc.t, t_g])
    n_o = len(lc.t)
    periods = hypS.period * np.exp(np.linspace(np.log(GRID_LO), np.log(GRID_HI), N_GRID))
    pg = Periodogram(lc.t, lc.e, lc.survey, lc.n_surveys, hypS, hypD, periods)
    dz_real = pg.dlnZ(lc.y[:, None])[N_GRID // 2, 0]       # at the published period
    acc_x, acc_dz, acc_P = [], [], []
    pool = []   # (dz, P, x) of all trials, kept for the relaxed fallback (top only)
    n_trials = 0
    while len(acc_x) < k and n_trials < N_TRIALS_MAX:
        if psd == "drw":
            X = M.simulate_drw(t_all, hypD.sigma, hypD.tau, rng, BATCH)
        else:
            X = simulate_dho(t_all, hypD.sigma, hypD.tau, rng, BATCH)
        Yo = (np.median(lc.y) + X[:, :n_o] + lc.e * rng.standard_normal((BATCH, n_o))
              + hypD.jitter * rng.standard_normal((BATCH, n_o))).T
        dz = pg.dlnZ(Yo)
        ib = np.argmax(dz, axis=0)
        best = dz[ib, np.arange(BATCH)]
        for b in np.where(best >= dz_real)[0]:
            acc_x.append((X[b], Yo[:, b].copy())); acc_dz.append(best[b]); acc_P.append(periods[ib[b]])
        top = np.argsort(best)[-k:]
        pool += [(best[b], periods[ib[b]], (X[b], Yo[:, b].copy())) for b in top]
        pool = sorted(pool, key=lambda z: -z[0])[:k]
        n_trials += BATCH
    relaxed = len(acc_x) < k
    if relaxed:
        acc_dz = [z[0] for z in pool]; acc_P = [z[1] for z in pool]; acc_x = [z[2] for z in pool]
    acc_rate = (len(acc_x) if not relaxed else 0) / n_trials
    rows = []
    for (x, yo), dzv, Pf in zip(acc_x[:k], acc_dz[:k], acc_P[:k]):
        y_o = yo                 # the same noisy realisation that passed the selection
        fake = M.LightCurve(lc.t, y_o, lc.e, lc.survey)
        hD = M.fit_hyper(fake, "D", n_starts=2, rng=rng)
        hS = M.fit_hyper(fake, "S", period=Pf, n_starts=2, rng=rng)
        s = np.exp(C.SCALE_PRIOR_SIGMA_LN * rng.standard_normal())
        j = abs(C.JITTER_GAIA_PRIOR_SIGMA * rng.standard_normal())
        c = OFFSET_SIM_SIGMA * rng.standard_normal()
        y_g = np.median(lc.y) + c + s * x[n_o:] + np.sqrt(e_g**2 + j**2) * rng.standard_normal(len(t_g))
        sc = S.CandidateScorer(fake, hS, hD, t_g, e_g)
        lnb = float(sc.lnB(y_g[None])[0])
        nul = sc.null_lnB(N_COND, rng)
        p = S.null_pvalue(lnb, nul)
        rows.append(dict(template=row["cand_id"], sample=row["sample"], P_template=hypS.period,
                         P_fake=Pf, dlnZ_fake=dzv, dlnZ_real=dz_real, lnB=lnb, p=p,
                         sigma_s_fit=hS.sigma_s, tau_fit=hD.tau, relaxed=relaxed,
                         acc_rate=acc_rate, n_trials=n_trials, psd=psd))
    return rows


# ---------------------------------------------------------------------------
def _template_job(args):
    cid, k, seed, psd = args
    import frozen  # noqa: WPS433  (local module: loads frozen fits + mock cadence)
    rng = np.random.default_rng(seed)
    row, lc, hS, hD, t_g, e_g = frozen.load_template(cid)
    return run_template(row, lc, hS, hD, t_g, e_g, k, rng, psd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=50)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--templates", default="")
    ap.add_argument("--psd", default="drw", choices=["drw", "dho"])
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    import frozen
    ids = a.templates.split(",") if a.templates else frozen.template_ids()
    out = Path(a.out) if a.out else C.RESULTS / f"null_sim_{a.psd}.csv"
    done = set(pd.read_csv(out)["template"]) if out.exists() else set()
    jobs = [(cid, a.k, C.RANDOM_SEED + 7919 * i + (0 if a.psd == "drw" else 1), a.psd)
            for i, cid in enumerate(ids) if cid not in done]
    print(f"{len(jobs)} templates to run, k={a.k}, workers={a.workers}", flush=True)
    t0 = time.time()
    if a.workers > 1:
        from multiprocessing import get_context
        with get_context("spawn").Pool(a.workers) as pool:
            for rows in pool.imap_unordered(_template_job, jobs):
                _append(out, rows, t0)
    else:
        for j in jobs:
            _append(out, _template_job(j), t0)


def _append(out, rows, t0):
    df = pd.DataFrame(rows)
    df.to_csv(out, mode="a", header=not out.exists(), index=False)
    print(f"{df.template.iloc[0]}: n={len(df)} relaxed={df.relaxed.iloc[0]} "
          f"acc={df.acc_rate.iloc[0]:.4f} frac(p<0.05)={np.mean(df.p < 0.05):.3f} "
          f"frac(lnB>3)={np.mean(df.lnB > 3):.3f}  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    main()
