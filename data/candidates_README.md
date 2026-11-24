# data/candidates.csv — published periodic-quasar candidates

Built by `src/build_candidates.py` from saved copies of the published tables, then checked
by `src/verify_positions.py`. Built 2026-11-24. **No ephemeris value is invented.** Values a
paper does not publish are blank and flagged in `flags`.

| Sample | Rows | Final sample | Source of rows | Period error | Amplitude | Phase / T0 |
|---|---|---|---|---|---|---|
| Graham2015 (`G15_*`) | 111 | 111 | VizieR J/MNRAS/453/1562 table2 | not published | not published | not published |
| Charisi2016 (`C16_*`) | 50 | 33 (bold in Table 1) | arXiv:1604.01020v2 LaTeX, Tables 1–2 | not published | derived from "variability fraction" | not published |
| Foustoul2025 (`F25_*`) | 257 | 36 (Table B.1) | B.1: arXiv PDF text; C.1: VizieR J/A+A/699/A55 | B.1: 3σ/3; C.1: published | B.1: ZTF g semi-amplitude | not published |
| Liu2019 (`L19_*`) | 26 | 1 (PSO J185.8689+46.9752) | arXiv PDF text, Tables 4, 5, 7 | published | PS1 r A0 | not published |
| Chen2024 (`C24_*`) | 3 | 3 | text, Sect. 3 | — | — | — (periods not in an accessible table) |

Verbatim copies of the parsed paper tables are in `data/paper_tables/`. arXiv v1 dates were
checked on arXiv; DOIs on Crossref.

## Columns
`cand_id, sample, name, ra_deg, dec_deg, coord_source, z, mag, mag_band, period_d (observed
frame), period_err_d, period_err_note, period_frame, period_band, ref_epoch_mjd, phase
(always blank: never published), amp_mag, amp_definition, in_final_sample, group
(forecast/gap/exploratory), arxiv_v1_date, source_bibcode, source_doi, source_arxiv,
source_table, foustoul2025_status, elbadry2026_tested, flags`, plus `charisi_pval_composite`,
`foustoul_Ar` and the position check `xm_catalog, xm_name, xm_sep_arcsec, xm_ra, xm_dec,
xm_z` (nearest SDSS DR16Q, otherwise Milliquas v8, within 5″; 25″ for truncated identifiers).

## Corrections and flags worth knowing
- Charisi+2016 Table 1 has three coordinate typos (J005453.30−003258.3 Dec "−0.0.5495320";
  J150900.70+175114.3 Dec "117.85"; J155308.65+501436.5 Dec "550.24"). Coordinates were taken
  from the IAU names (flag `table1_coord_typo`).
- Charisi amplitude: A = (5/ln 10)·asinh(VF/100), the inverse of the Table 2 note's formula
  (flag `amp_derived_from_VF`; peak-to-peak).
- Foustoul B.1 identifiers are truncated. 26 match Graham+2015 objects (as the paper states).
  The other 10 were matched to DR16Q/Milliquas (redshifts agree to ±0.01). Their positions are
  the catalogue positions (`coord_source` says so).
- `foustoul2025_status` for Graham rows: 26 "consistent with ZTF" (in B.1); HS 1630+2355
  "rejected" (named in the text); 84 "verdict not published" (Foustoul do not list the 58
  red-noise / 5 unclear objects).
- `elbadry2026_tested`: no candidate lies within 2″ of the 181 Gaia DR3 candidates tested by
  El-Badry, Hogg & Rix (2026) (list from IvS-KULeuven/GASP commit 6085247).
- No quasar-catalogue match within 5″ for 6 rows: HS 0423+0658 and 5 low-z Foustoul C.1
  objects (z = 0.20–0.47; probably galaxy-catalogue AGN). They are kept and flagged.
- z disagreements > 0.05 with DR16Q/Milliquas are flagged (`z_mismatch`); the paper value is
  kept.
