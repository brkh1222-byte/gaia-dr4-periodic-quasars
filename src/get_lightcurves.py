"""Download pre-Gaia (CRTS, PTF) and post-Gaia (ZTF) light curves for all candidates (T1.2).

Only ground-based survey data are downloaded here. NO Gaia photometry of any kind.

  CRTS  Catalina Surveys Data Release cone search (CSS 'photcat'), radius 0.1 arcmin
        http://nunuku.caltech.edu/cgi-bin/getcssconedb_release_imgN.cgi
  PTF   IRSA TAP table ptf_lightcurves (PTF + iPTF, Mould-R fid=2 and g fid=1), radius 2"
  ZTF   IRSA light-curve API (current public data release), radius 1.5"

Raw per-epoch files are written to data/lightcurves/<survey>/<cand_id>.csv together with a
log (data/lightcurves/<survey>_log.csv). The selection of the matching object ID and the
quality cuts are applied later, in src/lightcurves.py, so the raw files stay untouched.

Usage:  .venv/bin/python src/get_lightcurves.py crts|ptf|ztf [--samples Graham2015,...]
"""
from __future__ import annotations

import argparse
import io
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402

UA = {"User-Agent": "periodic-quasar-research-script/1.0"}
CRTS_URL = "http://nunuku.caltech.edu/cgi-bin/getcssconedb_release_imgN.cgi"
ZTF_URL = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
IRSA_TAP = "https://irsa.ipac.caltech.edu/TAP"


def _retry(fn, tries=4, wait=10):
    for k in range(tries):
        try:
            return fn()
        except Exception as e:  # network errors, timeouts, parse errors
            err = e
            time.sleep(wait * (k + 1))
    raise err


def get_crts(ra, dec):
    def call():
        data = {"RADec": f"{ra:.6f} {dec:+.6f}", "Rad": "0.1", "IMG": "nun", "DB": "photcat",
                "OUT": "csv", "SHORT": "short", ".submit": "Submit"}
        files = {k: (None, v) for k, v in data.items()}   # multipart/form-data
        r = requests.post(CRTS_URL, files=files, headers=UA, timeout=180)
        r.raise_for_status()
        m = re.search(r"(http://nunuku\.caltech\.edu/DataRelease/upload/[^\"'<> ]+\.csv)", r.text)
        if m is None:
            if re.search(r"no (objects|sources|data)|not found|0 objects", r.text, re.I):
                return pd.DataFrame()
            raise RuntimeError("CRTS: no CSV link in response")
        c = requests.get(m.group(1), headers=UA, timeout=180)
        c.raise_for_status()
        if not c.text.strip():
            return pd.DataFrame()
        return pd.read_csv(io.StringIO(c.text))
    return _retry(call)


def get_ptf(ra, dec, service=None):
    import pyvo
    service = service or pyvo.dal.TAPService(IRSA_TAP)
    q = ("SELECT oid, obsmjd, mag_autocorr, magerr_auto, fid, ra, dec, goodflag, photcalflag, "
         "flags, qaflags, limitmag FROM ptf_lightcurves WHERE CONTAINS(POINT('ICRS',ra,dec), "
         f"CIRCLE('ICRS',{ra:.6f},{dec:+.6f},{2.0/3600:.6f}))=1")
    return _retry(lambda: service.run_sync(q).to_table().to_pandas())


def get_ztf(ra, dec):
    def call():
        params = {"POS": f"CIRCLE {ra:.6f} {dec:+.6f} {1.5/3600:.6f}", "BANDNAME": "g,r,i",
                  "FORMAT": "csv"}
        r = requests.get(ZTF_URL, params=params, headers=UA, timeout=300)
        r.raise_for_status()
        return pd.read_csv(io.StringIO(r.text)) if r.text.strip() else pd.DataFrame()
    return _retry(call)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("survey", choices=["crts", "ptf", "ztf"])
    ap.add_argument("--samples", default="")
    ap.add_argument("--wait-hours", type=float, default=0.0,
                    help="CRTS only: probe the server every 15 min for up to this long")
    a = ap.parse_args()
    if a.survey == "crts" and a.wait_hours > 0:
        t_end = time.time() + 3600 * a.wait_hours
        while True:
            try:
                if len(get_crts(196.3875, -10.5553)) > 0:   # PG 1302-102 probe
                    print("CRTS server responding", flush=True)
                    break
            except Exception as e:
                print("CRTS probe failed:", str(e)[:80], pd.Timestamp.now("UTC"), flush=True)
            if time.time() > t_end:
                sys.exit("CRTS server still down; giving up")
            time.sleep(900)
    cand = pd.read_csv(C.CANDIDATES_CSV)
    if a.samples:
        cand = cand[cand["sample"].isin(a.samples.split(","))]
    cand = cand[np.isfinite(cand["ra_deg"])]
    outdir = C.LC_DIR / a.survey
    outdir.mkdir(parents=True, exist_ok=True)
    logf = C.LC_DIR / f"{a.survey}_log.csv"
    log = pd.read_csv(logf) if logf.exists() else pd.DataFrame(columns=["cand_id"])
    done = set(log["cand_id"])
    svc = None
    if a.survey == "ptf":
        import pyvo
        svc = pyvo.dal.TAPService(IRSA_TAP)
    for k, (_, r) in enumerate(cand.iterrows()):
        cid = r["cand_id"]
        if cid in done:
            continue
        t0 = time.time()
        try:
            if a.survey == "crts":
                df = get_crts(r["ra_deg"], r["dec_deg"])
            elif a.survey == "ptf":
                df = get_ptf(r["ra_deg"], r["dec_deg"], svc)
            else:
                df = get_ztf(r["ra_deg"], r["dec_deg"])
            status = "ok"
        except Exception as e:
            df, status = pd.DataFrame(), f"error: {str(e)[:120]}"
        if len(df):
            df.to_csv(outdir / f"{cid}.csv", index=False)
        rec = dict(cand_id=cid, n_rows=len(df), status=status,
                   retrieved_utc=pd.Timestamp.now("UTC").isoformat(timespec="seconds"),
                   seconds=round(time.time() - t0, 1))
        if status == "ok":          # failed queries are retried on the next run
            log = pd.concat([log, pd.DataFrame([rec])], ignore_index=True)
            log.to_csv(logf, index=False)
        print(f"[{k+1}/{len(cand)}] {cid} {status} n={len(df)}", flush=True)
        time.sleep(1.0)


if __name__ == "__main__":
    main()
