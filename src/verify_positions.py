"""Check every candidate position against quasar catalogues (task T1.1 provenance check).

For each row of data/candidates.csv we query (VizieR, CfA mirror):
  * SDSS DR16Q (Lyke et al. 2020; VizieR VII/289/dr16q)
  * Milliquas v8 (Flesch 2023; VizieR VII/294/catalog) as fallback (covers non-SDSS sky)
within 5" (25" for identifiers truncated to whole seconds). Results are cached in
data/raw/position_check.csv and merged into candidates.csv as xm_* columns. Rows whose
coordinates came only from a truncated identifier get the matched quasar position (flagged).
A redshift disagreement > 0.05 between the paper and the catalogue is flagged.
Columns ra_use/dec_use hold the position used downstream (catalogue position when matched,
otherwise the published one).
"""
from __future__ import annotations

import io
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402

ASU = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv"
CACHE = C.RAW / "position_check.csv"
CATS = {"DR16Q": ("VII/289/dr16q", "_r,RAJ2000,DEJ2000,z,SDSS"),
        "Milliquas": ("VII/294/catalog", "_r,RAJ2000,DEJ2000,z,Name,Type,Rmag")}


def cone(cat, ra, dec, radius):
    src, cols = CATS[cat]
    params = {"-source": src, "-c": f"{ra:.6f} {dec:+.6f}", "-c.rs": radius,
              "-out": cols, "-sort": "_r", "-out.max": 5}
    for attempt in range(4):
        try:
            r = requests.get(ASU, params=params, timeout=60)
            r.raise_for_status()
            break
        except requests.RequestException:
            time.sleep(5 * (attempt + 1))
    else:
        raise RuntimeError(f"VizieR query failed for {ra},{dec}")
    lines = [ln for ln in r.text.splitlines() if ln and not ln.startswith("#")]
    if len(lines) <= 3:
        return None
    df = pd.read_csv(io.StringIO("\n".join([lines[0]] + lines[3:])), sep="\t")
    df.columns = [c.strip() for c in df.columns]
    row = df.iloc[0]
    name = str(row.get("SDSS", row.get("Name", ""))).strip()
    return dict(xm_catalog=cat, xm_name=name, xm_sep_arcsec=float(row["_r"]),
                xm_ra=float(row["RAJ2000"]), xm_dec=float(row["DEJ2000"]),
                xm_z=float(row["z"]) if str(row["z"]).strip() not in ("", "nan") else np.nan)


def main():
    cand = pd.read_csv(C.CANDIDATES_CSV)
    cache = pd.read_csv(CACHE) if CACHE.exists() else pd.DataFrame(columns=["cand_id"])
    done = set(cache["cand_id"])
    new = []
    for _, r in cand.iterrows():
        if r["cand_id"] in done or not np.isfinite(r["ra_deg"]):
            continue
        trunc = "truncated" in str(r["flags"])
        rad = 25 if trunc else 5
        res = cone("DR16Q", r["ra_deg"], r["dec_deg"], rad) or \
            cone("Milliquas", r["ra_deg"], r["dec_deg"], rad) or {"xm_catalog": "none"}
        res["cand_id"] = r["cand_id"]
        new.append(res)
        if len(new) % 25 == 0:
            print(len(new), "queried", flush=True)
            cache = pd.concat([cache, pd.DataFrame(new)], ignore_index=True)
            cache.to_csv(CACHE, index=False); new = []
        time.sleep(0.3)
    cache = pd.concat([cache, pd.DataFrame(new)], ignore_index=True)
    cache.to_csv(CACHE, index=False)

    cand = cand.drop(columns=[c for c in cand.columns if c.startswith("xm_")])
    cand = cand.merge(cache, on="cand_id", how="left")
    for i, r in cand.iterrows():
        fl = str(r["flags"]) if isinstance(r["flags"], str) else ""
        if r["xm_catalog"] == "none" or pd.isna(r["xm_catalog"]):
            fl += ";no_quasar_catalogue_match"
        else:
            if "coords_from_truncated_name" in fl:
                cand.at[i, "ra_deg"], cand.at[i, "dec_deg"] = r["xm_ra"], r["xm_dec"]
                cand.at[i, "coord_source"] = f"{r['xm_catalog']} {r['xm_name']} (matched to truncated id, {r['xm_sep_arcsec']:.1f}\")"
                if not np.isfinite(r["z"]):
                    pass
            elif r["xm_sep_arcsec"] > 2.0:
                fl += f";catalogue_match_offset_{r['xm_sep_arcsec']:.1f}arcsec"
            if np.isfinite(r["z"]) and np.isfinite(r["xm_z"]) and abs(r["z"] - r["xm_z"]) > 0.05:
                fl += f";z_mismatch(paper={r['z']},{r['xm_catalog']}={r['xm_z']})"
        cand.at[i, "flags"] = fl.strip(";")
    # position used for all downloads and cross-matches: the catalogue position when matched
    m = cand["xm_catalog"].isin(["DR16Q", "Milliquas"])
    cand["ra_use"] = np.where(m, cand["xm_ra"], cand["ra_deg"])
    cand["dec_use"] = np.where(m, cand["xm_dec"], cand["dec_deg"])
    cand["pos_use_source"] = np.where(m, cand["xm_catalog"], "published")
    cand.to_csv(C.CANDIDATES_CSV, index=False, float_format="%.8g")
    print(cand["xm_catalog"].value_counts(dropna=False))
    print("offset>2\":", cand["flags"].str.contains("catalogue_match_offset").sum(),
          "z mismatch:", cand["flags"].str.contains("z_mismatch").sum(),
          "no match:", cand["flags"].str.contains("no_quasar").sum())


if __name__ == "__main__":
    main()
