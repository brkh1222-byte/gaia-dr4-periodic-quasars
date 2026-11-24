"""Assemble preregistration/ from the frozen results (task T2.4).

  * copies results/frozen/{frozen_params,frozen_predictions,power_mock}.csv and a snapshot of
    data/candidates.csv into preregistration/;
  * summarises results/null_sim_drw.csv (and _dho if present) into null_sim_summary.csv/json
    (r0, its 95% upper bound, r0_used = max(0.05, upper bound), false-survivor rate);
  * fills the {{...}} placeholders in preregistration/prereg.md (template kept as
    prereg_template.md);
  * writes CHECKSUMS.txt (SHA-256).

It never runs git: the researcher approves, then the commit is made by hand.
Usage: .venv/bin/python src/make_prereg.py --code-hash <commit> --power-rule 0.5
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402


def clopper_pearson_upper(k, n, cl=0.95):
    return 1.0 if k == n else float(beta.ppf(cl, k + 1, n - k))


def null_summary(path):
    d = pd.read_csv(path)
    n, k = len(d), int((d.p < C.P_SINGLE).sum())
    return dict(n_fakes=n, n_templates=int(d.template.nunique()),
                k_per_template=int(d.groupby("template").size().median()),
                r0=k / n, r0_upper95=clopper_pearson_upper(k, n),
                false_survivor_rate=float(np.mean(d.lnB > C.LNB_THRESHOLD)),
                frac_relaxed_templates=float(d.groupby("template").relaxed.first().mean()),
                median_lnB=float(d.lnB.median()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--code-hash", required=True)
    ap.add_argument("--power-rule", type=float, required=True)
    a = ap.parse_args()
    P, R = C.PREREG, C.RESULTS
    fr = R / "frozen"
    for f in ("frozen_params.csv", "frozen_predictions.csv", "power_mock.csv"):
        shutil.copy(fr / f, P / f)
    shutil.copy(C.CANDIDATES_CSV, P / "candidates_snapshot.csv")
    shutil.copy(C.RAW / "dr3_mean_g.csv", P / "dr3_mean_g_snapshot.csv")

    summ = {"drw": null_summary(R / "null_sim_drw.csv")}
    if (R / "null_sim_dho.csv").exists():
        summ["dho"] = null_summary(R / "null_sim_dho.csv")
    r0_used = max(C.P_SINGLE, summ["drw"]["r0_upper95"])
    summ["r0_used"] = r0_used
    json.dump(summ, open(P / "null_sim_summary.json", "w"), indent=2)
    json.dump(summ, open(R / "null_sim_summary.json", "w"), indent=2)
    pd.DataFrame(summ["drw"], index=["drw"]).to_csv(P / "null_sim_summary.csv")

    tmpl = P / "prereg_template.md"
    if not tmpl.exists():
        shutil.copy(P / "prereg.md", tmpl)
    text = tmpl.read_text()
    pw = pd.read_csv(fr / "power_mock.csv")
    fp = pd.read_csv(fr / "frozen_params.csv")
    ds = pd.read_csv(fr / "discovery_summary.csv")
    rows = []
    for s in ("Graham2015", "Charisi2016", "Foustoul2025"):
        tot = ds[ds["sample"] == s]
        sc = tot[tot.status == "scored"]
        p = pw[pw["sample"] == s]
        rows.append(f"| {s} | {len(tot)} | {len(sc)} | {int((p.power_p05 >= a.power_rule).sum())} | "
                    f"{p.power_p05.median():.2f} | {p.power_lnB3.median():.2f} |")
    table = ("| Sample | Rows considered | Scored (G<20.7, data OK) | Testable on mock cadence "
             f"(power ≥ {a.power_rule}) | Median power (p<0.05) | Median P(lnB>3 given S) |\n"
             "|---|---|---|---|---|---|\n" + "\n".join(rows))
    d = summ["drw"]
    repl = {
        "{{CODE_HASH}}": a.code_hash,
        "{{POWER_RULE}}": f"{a.power_rule}",
        "{{R0}}": f"{d['r0']:.4f}",
        "{{R0_USED}}": f"{r0_used:.4f}",
        "{{K}}": str(d["k_per_template"]),
        "{{NFAKE}}": str(d["n_fakes"]),
        "{{FALSE_SURV}}": f"{d['false_survivor_rate']:.4f}",
        "{{R0_DHO}}": f"{summ['dho']['r0']:.4f}" if "dho" in summ else "(not run)",
        "{{POWER_TABLE}}": table,
    }
    for k, v in repl.items():
        text = text.replace(k, v)
    (P / "prereg.md").write_text(text)
    # the testability threshold must already be frozen in the analysis code (no edits here)
    import re
    m = re.search(r"TESTABLE_POWER_P05 = ([0-9.]+)", (C.ROOT / "src" / "analyze_dr4.py").read_text())
    assert m and float(m.group(1)) == a.power_rule, "analyze_dr4.TESTABLE_POWER_P05 != --power-rule"

    files = sorted(f for f in P.iterdir() if f.is_file() and f.name != "CHECKSUMS.txt")
    with open(P / "CHECKSUMS.txt", "w") as fh:
        for f in files:
            fh.write(f"{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.name}\n")
    print(json.dumps(summ, indent=1))
    print("prereg.md filled; remaining placeholders:", [k for k in repl if k in text])


if __name__ == "__main__":
    main()
