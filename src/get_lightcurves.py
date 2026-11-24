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
import signal
import socket
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402

UA = {"User-Agent": "periodic-quasar-research-script/1.0"}
socket.setdefaulttimeout(120)
WATCHDOG_S = 600   # hard limit per object (covers DNS hangs that socket timeouts miss)


class _Watchdog(Exception):
    pass


def _alarm(signum, frame):
    raise _Watchdog("watchdog timeout")
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


CRTS_NAME_URL = "http://nunuku.caltech.edu/cgi-bin/getcssconedbid_release2.cgi"


def _crts_post(url, data, timeout=90):
    files = {k: (None, v) for k, v in data.items()}   # multipart/form-data
    r = requests.post(url, files=files, headers=UA, timeout=timeout)
    r.raise_for_status()
    m = re.search(r"(http://nunuku\.caltech\.edu/DataRelease/upload/[^\"'<> ]+\.csv)", r.text)
    if m is None:
        if re.search(r"There were 0 lines|no (objects|sources) found", r.text, re.I):
            return pd.DataFrame()          # genuine empty result
        raise RuntimeError("CRTS: no CSV link in response (service failure)")
    c = requests.get(m.group(1), headers=UA, timeout=timeout)
    c.raise_for_status()
    return pd.read_csv(io.StringIO(c.text)) if c.text.strip() else pd.DataFrame()


def get_crts(ra, dec, name=None, tries=2):
    """Cone search (0.1'); if the service fails, the name search (Sesame) as fallback.

    The CRTS position search has been intermittent since 2026-11-24; failures raise so the
    object is retried in a later pass.
    """
    last = None
    for _ in range(tries):
        try:
            return _crts_post(CRTS_URL, {"RADec": f"{ra:.6f} {dec:+.6f}", "Rad": "0.1",
                                         "IMG": "nun", "DB": "photcat", "OUT": "csv",
                                         "SHORT": "short", ".submit": "Submit"})
        except Exception as e:
            last = e
        if name:
            try:
                df = _crts_post(CRTS_NAME_URL, {"Name": name, "Rad": "0.1", "DB": "photcat",
                                                "OUT": "csv", "SHORT": "short", ".submit": "Submit"})
                if len(df):
                    df["query"] = f"name:{name}"
                return df
            except Exception as e:
                last = e
        time.sleep(20)
    raise last


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
    ap.add_argument("--shard", default="0/1", help="k/n: process every n-th candidate starting at k")
    ap.add_argument("--wait-hours", type=float, default=0.0,
                    help="keep making passes over failed objects for up to this long")
    a = ap.parse_args()
    cand = pd.read_csv(C.CANDIDATES_CSV)
    if a.samples:
        cand = cand[cand["sample"].isin(a.samples.split(","))]
    cand = cand[np.isfinite(cand["ra_use"])]
    k_sh, n_sh = (int(x) for x in a.shard.split("/"))
    cand = cand.iloc[k_sh::n_sh]
    outdir = C.LC_DIR / a.survey
    outdir.mkdir(parents=True, exist_ok=True)
    logf = C.LC_DIR / (f"{a.survey}_log.csv" if n_sh == 1 else f"{a.survey}_log_shard{k_sh}.csv")
    log = pd.read_csv(logf) if logf.exists() else pd.DataFrame(columns=["cand_id"])
    done = set()
    for lf in C.LC_DIR.glob(f"{a.survey}_log*.csv"):
        done |= set(pd.read_csv(lf)["cand_id"])
    svc = None
    if a.survey == "ptf":
        import pyvo
        svc = pyvo.dal.TAPService(IRSA_TAP)
    t_end = time.time() + 3600 * a.wait_hours
    passes = 0
    while True:
        passes += 1
        n_fail = _one_pass(a, cand, outdir, logf, log, done, svc)
        log = pd.read_csv(logf) if logf.exists() else log
        done |= set(log["cand_id"])
        print(f"pass {passes}: {n_fail} failures", flush=True)
        if n_fail == 0 or time.time() > t_end:
            break
        time.sleep(600)


def _one_pass(a, cand, outdir, logf, log, done, svc):
    n_fail = 0
    for k, (_, r) in enumerate(cand.iterrows()):
        cid = r["cand_id"]
        if cid in done:
            continue
        t0 = time.time()
        signal.signal(signal.SIGALRM, _alarm)
        signal.alarm(WATCHDOG_S)
        try:
            if a.survey == "crts":
                nm = r["name"] if not str(r["name"]).startswith("ZTF oid") else None
                df = get_crts(r["ra_use"], r["dec_use"], nm)
            elif a.survey == "ptf":
                df = get_ptf(r["ra_use"], r["dec_use"], svc)
            else:
                df = get_ztf(r["ra_use"], r["dec_use"])
            status = "ok"
        except (Exception, _Watchdog) as e:
            df, status = pd.DataFrame(), f"error: {str(e)[:120]}"
        finally:
            signal.alarm(0)
        if len(df):
            df.to_csv(outdir / f"{cid}.csv", index=False)
        rec = dict(cand_id=cid, n_rows=len(df), status=status,
                   retrieved_utc=pd.Timestamp.now("UTC").isoformat(timespec="seconds"),
                   seconds=round(time.time() - t0, 1))
        if status == "ok":          # failed queries are retried in the next pass
            log = pd.concat([log, pd.DataFrame([rec])], ignore_index=True)
            log.to_csv(logf, index=False)
        else:
            n_fail += 1
        print(f"[{k+1}/{len(cand)}] {cid} {status} n={len(df)}", flush=True)
        time.sleep(1.0)
    return n_fail


if __name__ == "__main__":
    main()
