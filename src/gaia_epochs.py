"""T3.2 (run ONLY after the pre-registration commit is public and DR4 is released).

Download Gaia DR4 epoch photometry for the cross-matched sources via astroquery DataLink in
batches (default 200 IDs), with retries and caching, into data/gaia_epochs/<cand_id>.csv.
Only the G band is kept (BP/RP are not used by the pre-registered analysis).

The data_release / retrieval_type strings are expected to be 'Gaia DR4' / 'EPOCH_PHOTOMETRY';
if the archive uses different strings on release day, pass them on the command line and
record them in an amendment before scoring.

Usage: .venv/bin/python src/gaia_epochs.py --i-confirm-prereg-is-public [--batch 200]
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
from gaia_xmatch import prereg_is_committed  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", type=int, default=200)
    ap.add_argument("--data-release", default="Gaia DR4")
    ap.add_argument("--retrieval-type", default="EPOCH_PHOTOMETRY")
    ap.add_argument("--i-confirm-prereg-is-public", action="store_true")
    a = ap.parse_args()
    if not (prereg_is_committed() and a.i_confirm_prereg_is_public):
        sys.exit("Refusing: the pre-registration must be committed and pushed first.")
    from astroquery.gaia import Gaia
    ids = pd.read_csv(C.DATA / "gaia_ids.csv").query("status == 'ok'")
    outdir = C.DATA / "gaia_epochs"; outdir.mkdir(parents=True, exist_ok=True)
    todo = ids[[not (outdir / f"{c}.csv").exists() for c in ids.cand_id]]
    sid2cid = dict(zip(todo.dr4_source_id.astype(int), todo.cand_id))
    sids = list(sid2cid)
    for i in range(0, len(sids), a.batch):
        chunk = sids[i:i + a.batch]
        for k in range(6):
            try:
                res = Gaia.load_data(ids=chunk, data_release=a.data_release,
                                     retrieval_type=a.retrieval_type, format="csv",
                                     data_structure="INDIVIDUAL")
                break
            except Exception as e:
                print("retry", k, str(e)[:100], flush=True); time.sleep(30 * (k + 1))
        else:
            print("batch failed", i, flush=True); continue
        for key, tables in res.items():
            for tab in tables:
                df = tab.to_pandas()
                if "source_id" not in df or len(df) == 0:
                    continue
                sid = int(df["source_id"].iloc[0])
                if "band" in df:
                    df = df[df["band"].astype(str).str.upper() == "G"]
                df.to_csv(outdir / f"{sid2cid.get(sid, sid)}.csv", index=False)
        print(f"batch {i // a.batch + 1}: {len(chunk)} ids", flush=True)


if __name__ == "__main__":
    main()
