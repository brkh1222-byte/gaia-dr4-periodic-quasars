"""Load, clean and bin survey light curves; define each candidate's discovery data and
Gaia test window (frozen rules; see preregistration/prereg.md).

Rules
-----
CRTS   keep the CRTS MasterID with most epochs whose median position is within 3" of the
       candidate; drop Magerr <= 0.
PTF    Mould-R (fid=2), goodflag=1 and photcalflag=1; oid with most R epochs.
ZTF    g band (zg), catflags == 0, magerr > 0; every oid with >= 20 epochs is kept and given
       its own mean level (one "survey" per oid).
All    epochs within 0.5 d of each other (same night) are averaged (inverse-variance
       weights); the error of a night is max(formal error, rms/sqrt(n)). Nights more than
       5 robust sigma (1.4826 MAD) from the light-curve median are dropped.

Discovery data and test window (MJD), per analysis group
  forecast (Graham+2015 CRTS; Charisi+2016 PTF/iPTF R):
      discovery = survey nights with t < cutoff, cutoff = arXiv v1 date of the paper;
      test window = [cutoff, end of DR4].
  gap (Foustoul+2025: CRTS + ZTF g):
      discovery = CRTS nights (any t) + ZTF g nights with t < arXiv v1 date;
      test window = (last CRTS night, first ZTF night) intersected with DR4.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import models as M  # noqa: E402

NIGHT_GAP_D = 0.5
MIN_ZTF_EPOCHS_PER_OID = 20
OUTLIER_SIGMA = 5.0


def mjd(datestr):
    return float(Time(datestr).mjd)


def bin_nights(t, y, e, gap=NIGHT_GAP_D):
    o = np.argsort(t)
    t, y, e = t[o], y[o], e[o]
    groups = np.concatenate([[0], np.cumsum(np.diff(t) > gap)])
    T, Y, E, N = [], [], [], []
    for g in np.unique(groups):
        m = groups == g
        w = 1.0 / e[m] ** 2
        yy = np.sum(w * y[m]) / w.sum()
        ee = np.sqrt(1.0 / w.sum())
        if m.sum() > 1:
            ee = max(ee, np.std(y[m], ddof=1) / np.sqrt(m.sum()))
        T.append(np.sum(w * t[m]) / w.sum()); Y.append(yy); E.append(ee); N.append(m.sum())
    return np.array(T), np.array(Y), np.array(E), np.array(N)


def clip_outliers(t, y, e, nsig=OUTLIER_SIGMA):
    med = np.median(y)
    mad = 1.4826 * np.median(np.abs(y - med))
    if mad <= 0:
        return t, y, e
    keep = np.abs(y - med) < nsig * np.sqrt(mad**2 + e**2)
    return t[keep], y[keep], e[keep]


def _nearest_id(df, idcol, racol, deccol, ra0, dec0, max_sep=3.0):
    best, nbest = None, -1
    for oid, g in df.groupby(idcol):
        dra = (np.median(g[racol]) - ra0) * np.cos(np.radians(dec0)) * 3600
        ddec = (np.median(g[deccol]) - dec0) * 3600
        if np.hypot(dra, ddec) <= max_sep and len(g) > nbest:
            best, nbest = oid, len(g)
    return best


def load_crts(cid, ra, dec):
    f = C.LC_DIR / "crts" / f"{cid}.csv"
    if not f.exists():
        return None
    df = pd.read_csv(f)
    df = df[(df["Magerr"] > 0) & (df["Mag"] > 0)]
    oid = _nearest_id(df, "MasterID", "RA", "Dec", ra, dec)
    if oid is None:
        return None
    g = df[df["MasterID"] == oid]
    t, y, e, _ = bin_nights(g["MJD"].values, g["Mag"].values, g["Magerr"].values)
    return clip_outliers(t, y, e)


def load_ptf(cid, ra, dec, fid=2):
    f = C.LC_DIR / "ptf" / f"{cid}.csv"
    if not f.exists():
        return None
    df = pd.read_csv(f)
    df = df[(df["fid"] == fid) & (df["goodflag"] == 1) & (df["photcalflag"] == 1)
            & (df["magerr_auto"] > 0)]
    oid = _nearest_id(df, "oid", "ra", "dec", ra, dec, max_sep=2.0)
    if oid is None:
        return None
    g = df[df["oid"] == oid]
    t, y, e, _ = bin_nights(g["obsmjd"].values, g["mag_autocorr"].values, g["magerr_auto"].values)
    return clip_outliers(t, y, e)


def load_ztf(cid, band="zg"):
    """Returns list of (t, y, e) per ZTF oid with enough epochs."""
    f = C.LC_DIR / "ztf" / f"{cid}.csv"
    if not f.exists():
        return []
    df = pd.read_csv(f)
    df = df[(df["filtercode"] == band) & (df["catflags"] == 0) & (df["magerr"] > 0)]
    out = []
    for oid, g in df.groupby("oid"):
        if len(g) < MIN_ZTF_EPOCHS_PER_OID:
            continue
        t, y, e, _ = bin_nights(g["mjd"].values, g["mag"].values, g["magerr"].values)
        out.append(clip_outliers(t, y, e))
    return out


def discovery_data(row):
    """(LightCurve, test_window(lo, hi), info dict) or (None, None, info) if unusable."""
    sample, cid = row["sample"], row["cand_id"]
    ra, dec = row["ra_deg"], row["dec_deg"]
    cutoff = mjd(row["arxiv_v1_date"])
    info = {"cand_id": cid}
    if sample == "Graham2015":
        lc = load_crts(cid, ra, dec)
        if lc is None:
            return None, None, {**info, "reason": "no CRTS light curve"}
        t, y, e = lc
        m = t < cutoff
        L = M.LightCurve(t[m], y[m], e[m])
        win = (cutoff, C.GAIA_DR4_END_MJD)
        info.update(survey="CRTS", n_disc=int(m.sum()), n_crts_after_cutoff_unused=int((~m).sum()))
    elif sample == "Charisi2016":
        lc = load_ptf(cid, ra, dec)
        if lc is None:
            return None, None, {**info, "reason": "no PTF R light curve"}
        t, y, e = lc
        m = t < cutoff
        L = M.LightCurve(t[m], y[m], e[m])
        win = (cutoff, C.GAIA_DR4_END_MJD)
        info.update(survey="PTF-R", n_disc=int(m.sum()))
    elif sample == "Foustoul2025":
        parts, labels = [], []
        crts = load_crts(cid, ra, dec)
        t_last_crts = -np.inf
        if crts is not None and len(crts[0]):
            parts.append(crts); t_last_crts = crts[0].max()
        ztf = [z for z in load_ztf(cid, "zg")]
        ztf = [(t[t < cutoff], y[t < cutoff], e[t < cutoff]) for t, y, e in ztf]
        ztf = [z for z in ztf if len(z[0]) >= MIN_ZTF_EPOCHS_PER_OID]
        if not ztf:
            return None, None, {**info, "reason": "no usable ZTF g light curve"}
        t_first_ztf = min(z[0].min() for z in ztf)
        # reference survey (index 0) = ZTF g oid with most epochs (Foustoul's reference band)
        ztf.sort(key=lambda z: -len(z[0]))
        allp = ztf + parts
        t = np.concatenate([p[0] for p in allp]); y = np.concatenate([p[1] for p in allp])
        e = np.concatenate([p[2] for p in allp])
        sv = np.concatenate([np.full(len(p[0]), k) for k, p in enumerate(allp)])
        L = M.LightCurve(t, y, e, sv)
        lo = max(t_last_crts, C.GAIA_DR4_START_MJD)
        win = (lo, min(t_first_ztf, C.GAIA_DR4_END_MJD))
        info.update(survey="ZTF-g" + ("+CRTS" if crts is not None else ""), n_disc=len(t),
                    n_ztf_oids=len(ztf))
    else:
        return None, None, {**info, "reason": "exploratory sample, not scored"}
    if L.t.size < 30:
        return None, None, {**info, "reason": f"fewer than 30 discovery nights ({L.t.size})"}
    if win[1] - win[0] < 180:
        return None, None, {**info, "reason": "test window shorter than 180 d"}
    info.update(t_disc_first=float(L.t.min()), t_disc_last=float(L.t.max()),
                win_lo=float(win[0]), win_hi=float(win[1]))
    return L, win, info
