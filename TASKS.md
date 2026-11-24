# TASKS: gaia-dr4-periodic-quasars

⏰ **Hard deadline: Phases 1–2 must be committed and pushed BEFORE Gaia DR4 is released on 2 Dec 2026.** Read PROPOSAL.md and CLAUDE.md first. Tick boxes as tasks complete.

## Phase 0: Setup
- [x] **T0.1** Create `requirements.txt` (numpy, scipy, pandas, astropy, astroquery, celerite2, matplotlib, requests, tqdm) and `src/config.py`. Done with `uv` (Python 3.11, `.venv`), versions pinned.

## Phase 1: Candidates and pre-Gaia data (before 2 Dec)
- [ ] **T1.1** Write `data/candidates.csv` with columns: id, RA, Dec, source paper, survey, period (d) ± error, reference epoch/phase, amplitude (mag), band, mean mag, and a flag for whether it is already tested by Foustoul 2025 / El-Badry 2026. Samples:
  - Graham+2015 CRTS (111);
  - Charisi+2016 PTF (33);
  - Foustoul+2025 (36 + 221);
  - Chen+2022/24 ZTF (3);
  - Liu+2019 PS1 (only those with G < 20.7).

  Use VizieR or the paper tables, and record the source of every row.
- [ ] **T1.2** Write `src/get_crts.py` and `src/get_ztf.py` to download CRTS DR2 and ZTF DR24 light curves for all candidates into `data/`.
- [ ] **T1.3** Write `notebooks/01_sample.ipynb` covering magnitude cuts and expected Gaia testability (G ≈ 20.7 limit).

## Phase 2: Models and PRE-REGISTRATION (before 2 Dec)
- [ ] **T2.1** Write `src/models.py`:
  - (a) sinusoid with the published period/phase/amplitude plus DRW noise;
  - (b) pure DRW (celerite2).

  Both are fitted only to pre-2014.5 data (CRTS/PTF). Each produces the predictive mean and covariance at arbitrary times in 2014–2020.
- [ ] **T2.2** Write `src/score.py`: the log predictive likelihood of new data under each model and the Bayes factor, with pre-specified free nuisance parameters only (Gaia G offset and scale relative to the original band).
- [ ] **T2.3** Write `src/simulate_null.py`: DRW-simulated fake "candidates" selected with a simple periodogram cut, to calibrate the Bayes factor distribution.
- [ ] **T2.4 PRE-REGISTRATION.** Write `preregistration/prereg.md` covering the hypotheses, sample cuts, frozen model parameters (`preregistration/frozen_params.csv`), scoring rule, the population test, and the analysis code commit hash. **The researcher reviews it; commit and push before 2 Dec 2026.** Never edit it afterwards; add dated amendments.

## Phase 3: After DR4 release (2 Dec onwards)
- [ ] **T3.1** Write `src/gaia_xmatch.py`: a TAP cross-match of candidates to DR4 `source_id` (1″ radius), saving `data/gaia_ids.csv`.
- [ ] **T3.2** Write `src/gaia_epochs.py`: batched `EPOCH_PHOTOMETRY` downloads (a few hundred IDs per call; retries; cache) into `data/gaia_epochs/`.
- [ ] **T3.3** Run the pre-registered scoring, then the population test vs the DRW null.
- [ ] **T3.4** Secondary (labelled exploratory): joint CRTS + Gaia + ZTF ~20-year fits for the survivors.
- [ ] **T3.5** Figures: Bayes factor distribution vs null; 3–4 example light curves with both predictions overlaid.

## Phase 4: Write-up
- [ ] **T4.1** Write `paper/notes.md` with bullet results, and fill in Results/Discussion in `paper/main.tex` (RNAAS template). **Claude drafts the paper; the researcher reviews and edits it; the paper carries an AI-use disclosure.** Authorship is undecided: keep the placeholders `[AUTHORS — TBD]` and `[AI-USE STATEMENT — TBD]`.
