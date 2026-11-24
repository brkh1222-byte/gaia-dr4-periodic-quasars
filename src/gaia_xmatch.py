"""T3.1 (run ONLY after the pre-registration commit is public and DR4 is released).

Cross-match the scored candidates to Gaia DR4 sources within 1 arcsec of the position used
downstream (ra_use/dec_use), writing data/gaia_ids.csv. Table name defaults to
'gaiadr4.gaia_source' and can be overridden on release day (recorded in an amendment).

Usage: .venv/bin/python src/gaia_xmatch.py [--table gaiadr4.gaia_source] [--i-confirm-prereg-is-public]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402

TAP = "https://gea.esac.esa.int/tap-server/tap"


def prereg_is_committed() -> bool:
    """Blindness guard: preregistration/prereg.md must be tracked and pushed."""
    try:
        tracked = subprocess.run(["git", "ls-files", "--error-unmatch", "preregistration/prereg.md"],
                                 cwd=C.ROOT, capture_output=True).returncode == 0
        ahead = subprocess.run(["git", "status", "-sb"], cwd=C.ROOT, capture_output=True,
                               text=True).stdout
        return tracked and "ahead" not in ahead.splitlines()[0]
    except Exception:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default="gaiadr4.gaia_source")
    ap.add_argument("--radius", type=float, default=1.0)
    ap.add_argument("--i-confirm-prereg-is-public", action="store_true")
    ap.add_argument("--test-dr3", action="store_true",
                    help="code test against gaiadr3.gaia_source (positions and mean G only)")
    a = ap.parse_args()
    if a.test_dr3:
        a.table = "gaiadr3.gaia_source"
    elif not (prereg_is_committed() and a.i_confirm_prereg_is_public):
        sys.exit("Refusing: the pre-registration must be committed and pushed first.")
    import pyvo
    svc = pyvo.dal.TAPService(TAP)
    fp = pd.read_csv(C.RESULTS / ("frozen_provisional" if a.test_dr3 else "frozen") / "frozen_params.csv")
    cand = pd.read_csv(C.CANDIDATES_CSV).set_index("cand_id")
    out = []
    for cid in fp.cand_id:
        r = cand.loc[cid]
        q = (f"SELECT source_id, ra, dec, phot_g_mean_mag, DISTANCE(POINT({r.ra_use:.7f},{r.dec_use:+.7f}),"
             f" POINT(ra,dec))*3600 AS sep FROM {a.table} WHERE 1=CONTAINS(POINT(ra,dec), "
             f"CIRCLE({r.ra_use:.7f},{r.dec_use:+.7f},{a.radius/3600:.7f})) ORDER BY sep")
        for k in range(6):
            try:
                t = svc.run_sync(q).to_table().to_pandas(); break
            except Exception:
                time.sleep(10 * (k + 1))
        else:
            out.append(dict(cand_id=cid, status="query_failed")); continue
        if len(t) == 0:
            out.append(dict(cand_id=cid, status="no_source")); continue
        out.append(dict(cand_id=cid, status="ok", dr4_source_id=int(t.source_id[0]),
                        sep_arcsec=float(t.sep[0]), n_within=len(t),
                        dr4_g_mean=float(t.phot_g_mean_mag[0])))
        time.sleep(0.2)
    pd.DataFrame(out).to_csv(C.DATA / ("gaia_ids_dr3test.csv" if a.test_dr3 else "gaia_ids.csv"), index=False)
    print(pd.DataFrame(out).status.value_counts())


if __name__ == "__main__":
    main()
