# STATUS — checkpoint 2026-11-24 (emergency shutdown)

## Done
- T0.1: `.venv` (uv, Python 3.11), pinned `requirements.txt`, `src/config.py`; CLAUDE.md/TASKS.md
  updated (Claude drafts the paper; researcher edits; AI-use disclosure; authorship TBD). Commit bf6134b.
- notes/novelty_check_2026-11-24.md: DR4 = 2 Dec 2026 ~12:00 CET (ESA); DR3 window ends 2017-05-28;
  no overlapping paper found (search not exhaustive; re-check ~30 Nov; follow up "Gaia Collaboration 2026" lead).
- src/models.py (T2.1 core): DRW and sinusoid+DRW as GPs with analytically marginalised basis
  (survey means, sinusoid a/b), period-node mixture weighted by discovery-data evidence, Gaia
  nuisance (offset analytic, scale and jitter by quadrature). tests/test_models.py: 12/12 pass,
  including brute-force checks of the predictive density.
- data/paper_tables/: verbatim text of Charisi 2016 T1/T2, Foustoul 2025 B.1, Liu 2019 T4/5/7/8, Chen 2024 excerpt.
- data/raw/vizier/: Graham 2015 table2 (111 rows) and Foustoul C.1 (221 rows), from the CfA VizieR mirror
  (CDS Strasbourg was unreachable): `https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/MNRAS/453/1562&-out.max=unlimited&-out.all`

## Key findings for the design
- No paper publishes a phase/T0. Graham and Charisi give no period error. "Published ephemeris"
  therefore = published period; phase and amplitude come from the discovery data (GP posterior), frozen.
- Charisi T1 coordinate typos: J005453 Dec "-0.0.5495320", J150900 Dec "117.85", J155308 Dec "550.24"
  -> fix from SDSS names, flag. The 33 = composite P<0.01 (30) + 3 with no composite value; confirm via bold in LaTeX.
- Charisi amplitude: A = (5/ln10)·asinh(VF/100) from "variability fraction" (derived, flag it).
- Foustoul B.1 names are truncated (need DR16Q/Milliquas match); errors are 3-sigma.
- Chen 2024: the 3 survivors' extended periods are not tabulated; full Table 1 is machine-readable only.
- Liu 2019: PS1 MDS light curves not public -> exploratory only.
- CRTS: form http://nunuku.caltech.edu/cgi-bin/getcssconedb_release_imgN.cgi works (POST RADec, Rad=0.1,
  DB=photcat, OUT=csv, SHORT=short; result CSV link in the response HTML). Data run to 2016, so cut at the cutoff.
- Proposed cutoff = arXiv v1 date of each ephemeris paper (Graham 2015-07-26; Charisi 2016-04-04).

## Next steps
1. Re-download arXiv sources (the earlier ones were truncated): `curl -sL -m 1500 https://arxiv.org/e-print/<id>`
   into data/raw/arxiv/<id>/ for 1604.01020, 2505.06656, 2206.11497, 1906.08315; check Charisi bold rows
   and look for Chen's full Table 1.
2. Write src/build_candidates.py -> data/candidates.csv (with provenance columns and flags), cross-match
   positions to SDSS DR16Q (VizieR VII/289 on the CfA mirror), Gaia DR3 mean G only (gaia_source; NO epoch photometry).
3. src/get_crts.py, src/get_ptf.py (IRSA TAP ptf_lightcurves), src/get_ztf.py (IRSA API).
4. src/score.py, src/simulate_null.py (fast null: Woodbury rank-2 sinusoid, fixed hyperparams; validate on a subset with refits),
   mock Gaia cadence; then notebooks/01_sample.ipynb, preregistration drafts (DO NOT git add), paper/main.tex.
Restart: `cd ~/Desktop/gaia-dr4-periodic-quasars && .venv/bin/python -m pytest -q tests`
