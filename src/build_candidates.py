"""Build data/candidates.csv from the published candidate tables (tasks T1.1).

Every row records where each value came from. No ephemeris value is invented:
values that a paper does not publish are left blank and flagged.

Sources (all saved under data/raw/vizier or data/paper_tables before parsing):
  Graham+2015   VizieR J/MNRAS/453/1562 table2 (111 rows)
  Charisi+2016  arXiv:1604.01020v2 LaTeX source, Tables 1-2 (50 rows; 33 bold = final sample)
  Foustoul+2025 Table B.1 (36, arXiv PDF text) and VizieR J/A+A/699/A55 tablec1 (221)
  Liu+2019      Tables 4, 5 and 7 (26 rows, arXiv PDF text)
  Chen+2024     the 3 candidates that keep periodicity with extended data (text, Sect. 3)

Usage:  .venv/bin/python src/build_candidates.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402

VZ = C.RAW / "vizier"
PT = C.DATA / "paper_tables"

# arXiv v1 submission dates (checked on the arXiv abstract pages, 2026-11-24)
ARXIV_V1 = {
    "Graham2015": ("1507.07603", "2015-07-27"),
    "Charisi2016": ("1604.01020", "2016-04-04"),
    "Liu2019": ("1906.08315", "2019-06-19"),
    "Chen2024": ("2206.11497", "2022-06-23"),
    "Foustoul2025": ("2505.06656", "2025-05-10"),
}
BIBCODE = {
    "Graham2015": "2015MNRAS.453.1562G",
    "Charisi2016": "2016MNRAS.463.2145C",
    "Liu2019": "2019ApJ...884...36L",
    "Chen2024": "2024MNRAS.527.12154C",
    "Foustoul2025": "2025A&A...699A..55F",
}
DOI = {
    "Graham2015": "10.1093/mnras/stv1726",
    "Charisi2016": "10.1093/mnras/stw1838",
    "Liu2019": "10.3847/1538-4357/ab40cb",
    "Chen2024": "10.1093/mnras/stad3981",
    "Foustoul2025": "10.1051/0004-6361/202452818",
}
# Analysis group (see preregistration): forecast = discovery data end before the Gaia
# test epochs; gap/backcast = Gaia epochs fall in a gap of the discovery data;
# exploratory = cannot be tested as pre-registered (no public discovery light curve or
# no tabulated period).
GROUP = {"Graham2015": "forecast", "Charisi2016": "forecast", "Foustoul2025": "gap",
         "Liu2019": "exploratory", "Chen2024": "exploratory"}

COLUMNS = [
    "cand_id", "sample", "name", "ra_deg", "dec_deg", "coord_source", "z",
    "mag", "mag_band", "period_d", "period_err_d", "period_err_note", "period_frame",
    "period_band", "ref_epoch_mjd", "phase", "amp_mag", "amp_definition",
    "in_final_sample", "group", "arxiv_v1_date", "source_bibcode", "source_doi",
    "source_arxiv", "source_table", "foustoul2025_status", "flags",
]


def read_vizier_tsv(path: Path) -> pd.DataFrame:
    lines = [ln.rstrip("\n") for ln in path.open() if not ln.startswith("#") and ln.strip()]
    header = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[3:]]  # skip units and dashes
    df = pd.DataFrame(rows, columns=header)
    return df.apply(lambda s: s.str.strip())


def coord_from_sdss_name(name: str):
    """Decimal RA/Dec from an IAU 'Jhhmmss.ss+ddmmss.s' designation (or None)."""
    m = re.search(r"J(\d{2})(\d{2})(\d{2}(?:\.\d+)?)([+-])(\d{2})(\d{2})(\d{2}(?:\.\d+)?)", name)
    if not m:
        return None
    h, mi, s, sg, d, dm, ds = m.groups()
    c = SkyCoord(f"{h}h{mi}m{s}s", f"{sg}{d}d{dm}m{ds}s", frame="icrs")
    return c.ra.deg, c.dec.deg


def base_row(sample):
    arx, v1 = ARXIV_V1[sample]
    return dict(sample=sample, arxiv_v1_date=v1, source_bibcode=BIBCODE[sample],
                source_doi=DOI[sample], source_arxiv=arx, group=GROUP[sample],
                period_frame="observed", ref_epoch_mjd=np.nan, phase=np.nan, flags="")


def add_flag(row, flag):
    row["flags"] = (row["flags"] + ";" + flag).strip(";")


# ---------------------------------------------------------------------------
def graham2015():
    df = read_vizier_tsv(VZ / "J_MNRAS_453_1562.tsv")
    out = []
    for _, r in df.iterrows():
        row = base_row("Graham2015")
        c = SkyCoord(r["RAJ2000"], r["DEJ2000"], unit=(u.hourangle, u.deg))
        row.update(name=r["Name"], ra_deg=c.ra.deg, dec_deg=c.dec.deg,
                   coord_source="VizieR J/MNRAS/453/1562 table2", z=float(r["z"]),
                   mag=float(r["Vmed"]), mag_band="CRTS V (median)",
                   period_d=float(r["Per"]), period_err_d=np.nan,
                   period_err_note="not published", period_band="CRTS unfiltered (V-calibrated)",
                   amp_mag=np.nan, amp_definition="not published", in_final_sample=True,
                   source_table="Table 2 (VizieR J/MNRAS/453/1562/table2)")
        add_flag(row, "no_period_error_published")
        add_flag(row, "no_phase_published")
        add_flag(row, "no_amplitude_published")
        out.append(row)
    return out


# ---------------------------------------------------------------------------
def _strip_tex(s):
    s = re.sub(r"\$\^\{[^}]*\}\$", "", s)          # footnote markers
    s = s.replace("\\textbf{", "").replace("}", "").replace("$", "")
    return s.strip()


def charisi2016():
    tex = (PT / "charisi2016_tables.tex").read_text()
    t1, t2 = tex.split("Table:BinaryProperties", 1) if "Table:BinaryProperties" in tex else (tex, "")
    # Table 1 rows
    rows1 = {}
    for ln in t1.splitlines():
        if "&" not in ln or ln.strip().startswith("Name"):
            continue
        cells = [c.strip() for c in ln.split("&")]
        if len(cells) < 5:
            continue
        raw_name = cells[0]
        bold = "textbf" in raw_name
        name = _strip_tex(raw_name)
        try:
            z = float(cells[3]); rmag = float(cells[4])
        except ValueError:
            continue
        rows1[name] = dict(bold=bold, ra_raw=cells[1], dec_raw=cells[2], z=z, rmag=rmag)
    # Table 2 rows: Name & P-value & Period & VarFrac & ...
    rows2 = {}
    for ln in t2.splitlines():
        if "&" not in ln:
            continue
        cells = [c.strip() for c in ln.split("&")]
        name = _strip_tex(cells[0])
        try:
            per = float(cells[2]); vf = float(cells[3])
        except (ValueError, IndexError):
            continue
        rows2[name] = dict(pval=_strip_tex(cells[1]), period=per, vf=vf)
    out = []
    for name, a in rows1.items():
        row = base_row("Charisi2016")
        b = rows2.get(name)
        if b is None:
            raise RuntimeError(f"Charisi: no Table 2 row for {name}")
        flags = []
        ra, dec = a["ra_raw"], a["dec_raw"]
        try:
            ra_f, dec_f = float(ra), float(dec)
            if not (-90 <= dec_f <= 90):
                raise ValueError
            src = "Charisi+2016 Table 1"
        except ValueError:
            ra_f, dec_f = np.nan, np.nan
        nc = coord_from_sdss_name(name)
        if nc is not None:
            sep = (SkyCoord(nc[0], nc[1], unit="deg").separation(
                SkyCoord(ra_f if np.isfinite(ra_f) else 0, dec_f if np.isfinite(dec_f) else 0,
                         unit="deg")).arcsec)
            if not np.isfinite(ra_f) or sep > 2.0:
                flags.append(f"table1_coord_typo(RA={ra},Dec={dec});coords_from_IAU_name")
                ra_f, dec_f = nc
                src = "IAU name (Table 1 value is a typo)"
        # amplitude from the 'variability fraction' VF = sinh(A ln10 / 5) (Table 2 note)
        amp = 5.0 / np.log(10.0) * np.arcsinh(b["vf"] / 100.0)
        row.update(name=name, ra_deg=ra_f, dec_deg=dec_f, coord_source=src, z=a["z"],
                   mag=a["rmag"], mag_band="R (average, HMQ)", period_d=b["period"],
                   period_err_d=np.nan, period_err_note="not published",
                   period_band="PTF Mould-R", amp_mag=round(amp, 4),
                   amp_definition=("peak-to-peak A of best-fit sinusoid, derived from Table 2 "
                                   "'variability fraction' VF via A=(5/ln10)asinh(VF)"),
                   in_final_sample=a["bold"],
                   source_table="Tables 1-2 (arXiv source); final 33 shown in bold")
        for f in flags + ["no_period_error_published", "no_phase_published", "amp_derived_from_VF"]:
            add_flag(row, f)
        if b["pval"].startswith("—") or b["pval"].startswith("-"):
            add_flag(row, "no_extended_data(composite P-value not given)")
        row["charisi_pval_composite"] = b["pval"]
        out.append(row)
    return out


# ---------------------------------------------------------------------------
B1_RE = re.compile(
    r"^(?P<id>\S+(?: \S+)?) (?P<P>\d+) ± (?P<eP>\d+) (?P<chi>[\d.]+) \((?P<dof>\d+)\) "
    r"(?P<Pgls>\d+) / (?P<pgls>[\d.]+) (?P<Ar>[\d.]+) ± (?P<eAr>[\d.]+) "
    r"(?P<Ag>[\d.]+) ± (?P<eAg>[\d.]+) (?P<rest>.*)$")


def foustoul2025(graham_rows):
    out = []
    gcoords = SkyCoord([r["ra_deg"] for r in graham_rows], [r["dec_deg"] for r in graham_rows],
                       unit="deg")
    gnames = [r["name"] for r in graham_rows]
    for ln in (PT / "foustoul2025_tableB1.txt").read_text().splitlines():
        m = B1_RE.match(ln.strip())
        if not m:
            continue
        d = m.groupdict()
        row = base_row("Foustoul2025")
        ident = d["id"]
        # resolve the identifier: Graham name, or truncated Jhhmmss+ddmmss
        flags, src = [], ""
        gmatch = [i for i, n in enumerate(gnames) if n.replace(" ", "") == ident.replace(" ", "")]
        ra = dec = np.nan
        if gmatch:
            ra, dec = graham_rows[gmatch[0]]["ra_deg"], graham_rows[gmatch[0]]["dec_deg"]
            src = f"Graham+2015 position of {gnames[gmatch[0]]}"
        else:
            mm = re.match(r"J(\d{2})(\d{2})(\d{2})([+-])(\d{2})(\d{2})(\d{2})$", ident)
            if mm:
                h, mi, s, sg, dd, dm, ds = mm.groups()
                c = SkyCoord(f"{h}h{mi}m{s}.5s", f"{sg}{dd}d{dm}m{ds}.5s")
                sep = c.separation(gcoords).arcsec
                if sep.min() < 20:
                    i = int(np.argmin(sep))
                    ra, dec = graham_rows[i]["ra_deg"], graham_rows[i]["dec_deg"]
                    src = f"Graham+2015 position of {gnames[i]} (truncated name, {sep.min():.1f}\")"
                else:
                    ra, dec = c.ra.deg, c.dec.deg
                    src = "truncated identifier (needs catalogue match)"
                    flags.append("coords_from_truncated_name")
            else:
                flags.append("identifier_unresolved")
        P, eP3 = float(d["P"]), float(d["eP"])
        rest = d["rest"].split()
        row.update(name=ident, ra_deg=ra, dec_deg=dec, coord_source=src, z=np.nan,
                   mag=np.nan, mag_band="", period_d=P, period_err_d=round(eP3 / 3.0, 2),
                   period_err_note="published as 3-sigma; divided by 3",
                   period_band="CRTS (intercalibrated to ZTF g) + ZTF g",
                   amp_mag=float(d["Ag"]),
                   amp_definition="ZTF g semi-amplitude A in y=A sin(2pi t/P+t0)+B (Table B.1 col 6)",
                   in_final_sample=True, source_table="Table B.1 (36 candidates)")
        row["foustoul_Ar"] = float(d["Ar"])
        # z is the second-to-last-but-two field when separation is present
        try:
            nums = [x for x in rest]
            # rest = [sep] z logM ± elogM Pr ± ePr R  (sep missing for some rows)
            if "±" in nums:
                k = nums.index("±")
                row["z"] = float(nums[k - 2]) if k >= 2 else float(nums[0])
                if k == 1:  # no separation and no logM: rest = z Pr ± ePr R
                    row["z"] = float(nums[0])
        except (ValueError, IndexError):
            flags.append("z_parse_failed")
        if src.startswith("Graham"):
            gi = gnames.index(src.split("position of ")[1].split(" (")[0])
            zg = graham_rows[gi]["z"]
            if not (np.isfinite(row["z"]) and abs(row["z"] - zg) < 0.02):
                flags.append(f"z_from_Graham(TableB1_parse={row['z']})")
            row["z"] = zg
        for f in flags + ["no_phase_published(t0 fitted, origin not stated)"]:
            add_flag(row, f)
        out.append(row)
    if len(out) != 36:
        raise RuntimeError(f"Foustoul B.1 parsed {len(out)} rows, expected 36")
    # Table C.1 (VizieR)
    df = read_vizier_tsv(VZ / "J_A_A_699_A55.tsv")
    for _, r in df.iterrows():
        row = base_row("Foustoul2025")
        row.update(name=f"ZTF oid {r['Name']}", ra_deg=float(r["RAJ2000"]),
                   dec_deg=float(r["DEJ2000"]), coord_source="VizieR J/A+A/699/A55 tablec1",
                   z=float(r["z"]), mag=np.nan, mag_band="", period_d=float(r["Per"]),
                   period_err_d=float(r["e_Per"]), period_err_note="published (1-sigma per VizieR 'Period error')",
                   period_band="ZTF g (+CRTS if CRTSobs=Yes)", amp_mag=np.nan,
                   amp_definition="not published", in_final_sample=False,
                   source_table="Table C.1 weaker candidates (VizieR tablec1)")
        add_flag(row, "weaker_candidate_tableC1")
        add_flag(row, f"CRTSobs={r['CRTSobs']}")
        add_flag(row, "no_phase_published")
        out.append(row)
    return out


# ---------------------------------------------------------------------------
def liu2019():
    txt = (PT / "liu2019_tables.txt").read_text().replace("–", "-")
    t4 = dict()
    for m in re.finditer(r"PSO J([\d.]+)([+-][\d.]+) (\d+)±(\d+) \(", txt):
        t4[(m.group(1), m.group(2))] = (float(m.group(3)), float(m.group(4)))
    t5 = dict()
    for m in re.finditer(r"PSO J([\d.]+)([+-][\d.]+) \(([\d., ]+)\) \(([\d., ]+)\)", txt):
        mags = [float(x) for x in m.group(3).split(",")]
        amps = [float(x) for x in m.group(4).split(",")]
        t5[(float(m.group(1)), float(m.group(2)))] = (mags, amps)
    t7 = dict()
    for ln in txt.splitlines():
        m = re.match(r"PSO J([\d.]+)([+-][\d.]+) (SDSS|GS\d+[AB]|DCT\d+Q\d) ", ln)
        if m:
            toks = ln.split()
            # ... log(MBH) z Prest a a  -> z is 4th from the end
            try:
                t7[(float(m.group(1)), float(m.group(2)))] = float(toks[-4])
            except ValueError:
                pass
    out = []
    for (ra_s, dec_s), (P, eP) in t4.items():
        row = base_row("Liu2019")
        def fuzzy(tab):
            for (ra_k, de_k), v in tab.items():
                if abs(ra_k - float(ra_s)) < 2e-4 and abs(de_k - float(dec_s)) < 2e-4:
                    return v
            return None
        mags, amps = fuzzy(t5) or ([np.nan] * 4, [np.nan] * 4)
        zz = fuzzy(t7)
        row.update(name=f"PSO J{ra_s}{dec_s}", ra_deg=float(ra_s), dec_deg=float(dec_s),
                   coord_source="PS1 designation (Liu+2019 Table 4)", z=zz if zz is not None else np.nan,
                   mag=mags[1], mag_band="PS1 r (mean, Table 5)", period_d=P, period_err_d=eP,
                   period_err_note="published (Table 4, P_LS ± ΔP)", period_band="PS1 griz (filter-averaged)",
                   amp_mag=amps[1], amp_definition="PS1 r A0 of sinusoidal fit (Table 5; semi-amplitude implied, not stated)",
                   in_final_sample=(ra_s == "185.8689"),
                   source_table="Tables 4, 5, 7")
        add_flag(row, "no_phase_published")
        add_flag(row, "discovery_lc_not_public(PS1 MDS)")
        if fuzzy(t5) is None:
            add_flag(row, "table5_row_missing")
        if zz is None:
            add_flag(row, "table7_row_missing")
        out.append(row)
    if len(out) != 26:
        raise RuntimeError(f"Liu: parsed {len(out)} rows, expected 26")
    return out


# ---------------------------------------------------------------------------
def chen2024():
    out = []
    for name in ["SDSS J083111.23+123107.8", "SDSS J155304.40+240255.1",
                 "SDSS J135830.46-012255.7"]:
        row = base_row("Chen2024")
        ra, dec = coord_from_sdss_name(name)
        row.update(name=name, ra_deg=ra, dec_deg=dec, coord_source="IAU name", z=np.nan,
                   mag=np.nan, mag_band="", period_d=np.nan, period_err_d=np.nan,
                   period_err_note="", period_band="ZTF r (+CRTS, PTF, PS1, ATLAS extension)",
                   amp_mag=np.nan, amp_definition="not published", in_final_sample=True,
                   source_table="Sect. 3 text / Fig. 11 (3 candidates keeping periodicity)")
        add_flag(row, "period_not_in_accessible_table(ZTF-only P in machine-readable Table 1; extended P only in Fig. 11)")
        add_flag(row, "no_phase_published")
        out.append(row)
    return out


# ---------------------------------------------------------------------------
def foustoul_status(rows):
    """Flag Graham rows with Foustoul+2025 verdicts that can be identified."""
    fb = [r for r in rows if r["sample"] == "Foustoul2025" and r["in_final_sample"]]
    fc = SkyCoord([r["ra_deg"] for r in fb], [r["dec_deg"] for r in fb], unit="deg")
    for r in rows:
        if r["sample"] != "Graham2015":
            r["foustoul2025_status"] = ""
            continue
        sep = SkyCoord(r["ra_deg"], r["dec_deg"], unit="deg").separation(fc).arcsec
        if sep.min() < 3:
            r["foustoul2025_status"] = "consistent with ZTF (in Table B.1)"
        elif r["name"].replace(" ", "") == "HS1630+2355":
            r["foustoul2025_status"] = "rejected (named in text)"
        else:
            r["foustoul2025_status"] = "verdict not published (58 red-noise/5 unclear/22 no ZTF, unnamed)"


def main():
    g = graham2015()
    rows = g + charisi2016() + foustoul2025(g) + liu2019() + chen2024()
    foustoul_status(rows)
    df = pd.DataFrame(rows)
    df.insert(0, "cand_id", [f"{s[:2].upper()}{i:03d}" for s, i in
                             zip(df["sample"], df.groupby("sample").cumcount() + 1)])
    extra = [c for c in df.columns if c not in COLUMNS]
    df = df[COLUMNS + extra]
    df.to_csv(C.CANDIDATES_CSV, index=False, float_format="%.7g")
    print(df.groupby("sample").agg(n=("cand_id", "size"), final=("in_final_sample", "sum")))
    print("flags with coord problems:", df["flags"].str.contains("coord").sum())
    return df


if __name__ == "__main__":
    main()
