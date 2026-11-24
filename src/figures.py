"""Figures (T3.5). Fig. 1: ln B of testable candidates vs stacked null and injection
distributions (left); posterior of the fraction f of true periodic sources (right).

Usage: .venv/bin/python src/figures.py results/<tag>   (writes results/<tag>/fig1.pdf/.png)

The null/injection draws are regenerated with the same seeds as analyze_dr4.py, so the figure
uses exactly the distributions behind the p-values.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import analyze_dr4 as A  # noqa: E402
import frozen as F  # noqa: E402
import gaia as G  # noqa: E402
import lightcurves as LCm  # noqa: E402
import score as S  # noqa: E402


def regenerate_draws(tag_dir, epoch_dir, n_sim=1000):
    sc = pd.read_csv(tag_dir / "scores.csv")
    T = sc[(sc.status == "scored") & (sc.testable == True)]  # noqa: E712
    fp = pd.read_csv(C.RESULTS / "frozen" / "frozen_params.csv").set_index("cand_id")
    cand = pd.read_csv(C.CANDIDATES_CSV).merge(pd.read_csv(C.RAW / "dr3_mean_g.csv"),
                                                 on="cand_id").set_index("cand_id")
    nulls, injs = [], []
    for cid in T.cand_id:
        row = cand.loc[cid].copy(); row["cand_id"] = cid
        lc, win, _ = LCm.discovery_data(row)
        t_g, y_g, e_g = G.clean_dr4_epochs(pd.read_csv(epoch_dir / f"{cid}.csv"), win)
        hS, hD = A.load_frozen(fp.loc[cid])
        s = S.CandidateScorer(lc, hS, hD, t_g, e_g)
        rng = np.random.default_rng(C.RANDOM_SEED + F.cid_seed(cid) % 100000)
        nulls.append(s.null_lnB(n_sim, rng)); injs.append(s.inj_lnB(n_sim, rng))
    return T, nulls, injs


def fig1(tag_dir, epoch_dir, title=""):
    T, nulls, injs = regenerate_draws(tag_dir, epoch_dir)
    post = S.fraction_posterior(T.lnB.values, nulls, injs)
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
    lo, hi = -30, 40
    bins = np.linspace(lo, hi, 50)
    clip = lambda x: np.clip(x, lo, hi)  # noqa: E731
    w = lambda arrs: np.concatenate([np.full(len(a), len(T) / (len(a) * len(arrs))) for a in arrs])  # noqa: E731
    ax[0].hist(clip(np.concatenate(nulls)), bins, weights=w(nulls), histtype="step", color="0.4",
               lw=1.5, label="red-noise null (stacked)")
    ax[0].hist(clip(np.concatenate(injs)), bins, weights=w(injs), histtype="step", color="C1",
               lw=1.5, ls="--", label="if ephemerides true")
    ax[0].hist(clip(T.lnB.values), bins, color="C0", alpha=0.6, label=f"candidates (n={len(T)})")
    ax[0].axvline(C.LNB_THRESHOLD, color="k", lw=0.8, ls=":")
    ax[0].set_xlabel(r"predictive $\ln B$ (sinusoid vs. DRW)"); ax[0].set_ylabel("candidates per bin")
    ax[0].legend(fontsize=7, frameon=False)
    ax[1].plot(post["f"], post["post"], color="C0")
    ax[1].axvline(post["upper95"], color="k", ls=":", lw=0.8)
    ax[1].text(post["upper95"], ax[1].get_ylim()[1] * 0.9, f"  95%: {post['upper95']:.2f}", fontsize=7)
    ax[1].set_xlabel("fraction of true periodic sources $f$"); ax[1].set_ylabel("posterior density")
    if title:
        fig.suptitle(title, fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(tag_dir / f"fig1.{ext}", dpi=200)
    return post


if __name__ == "__main__":
    tag = Path(sys.argv[1])
    kind = tag.name.replace("dress_rehearsal_", "")
    epoch_dir = C.DATA / (f"mock_epochs_{kind}" if tag.name.startswith("dress") else "gaia_epochs")
    p = fig1(tag, epoch_dir, title=tag.name)
    print(json.dumps({k: v for k, v in p.items() if k in ("median", "upper95")}))
