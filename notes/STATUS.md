# STATUS — 2026-11-24 ~19:15 UTC — STOP POINT REACHED (awaiting approval of preregistration/)

## Where this project is (28 Nov 2026, plain language)

**About 60% done overall. Everything that can be done before the Gaia release is done or running;
the one thing blocking is your approval of the pre-registration.**

- **Finished:** the candidate list (447 objects from 5 papers, every value sourced and checked);
  the light curves we could get (PTF, ZTF, part of CRTS); the two models and the scoring code
  (24 automatic tests pass); the red-noise calibration (13,400 simulated fake candidates); a full
  practice run on fake Gaia data, which gave the right answer in both directions; the scripts
  that will download the real Gaia data on release day (tested on old Gaia data for stars that
  are not candidates); a paper draft with Introduction, Data and Methods.
- **Running:** retries of the Catalina (CRTS) light-curve server, which has been failing for most
  objects since 24 Nov. Only 38 of the 111 Graham candidates have Catalina data.
- **Next (before 2 Dec):** refresh the frozen predictions to include 12 more Foustoul objects;
  repeat the "has anyone done this already?" search; final practice run; polish the paper.
- **Needs you:** (1) say OK to the pre-registration (one page: notes/prereg_summary.md). It must
  be public before 2 Dec 12:00 CET or the test is no longer blind. (2) Confirm three choices:
  the 50%-power testability rule; using PTF data as a separate secondary group for Graham objects
  without Catalina data; merging the work branch into `main` before the pre-registration commit.
- **After 2 Dec:** download Gaia DR4 for the candidates, run the frozen analysis (about a day),
  fill in Results and Discussion, hand you the full draft.

## Decisions
- **Researcher, option A:** freeze each published period exactly as printed; phase and
  amplitude from discovery data before each paper's arXiv v1 date (GP posterior); frozen in
  preregistration/frozen_params.csv before DR4. State in prereg and paper that no source
  published a phase/T0 and that Graham/Charisi gave no period errors.
- Conditional null: per-candidate p-value from model-D predictive draws given the real
  discovery data (selection is then ignorable). End-to-end selected-DRW fake simulation
  (simulate_null.py) gives r0 for the H1 binomial rate, before release.

## Done (pushed)
- T0.1, T1.1 (data/candidates.csv: 447 rows; positions verified, 441 matched).
- Code: models, score, lightcurves, gaia, frozen, simulate_null, analyze_dr4,
  make_prereg, build_candidates, verify_positions, get_lightcurves; 21 tests pass.
- Paper skeleton paper/main.tex + verified refs.bib (no TeX installed locally).
- PTF light curves for Charisi (50/50).

## Running (background) — restart commands
- CRTS (intermittent service; empty answers now retried):
  `nice .venv/bin/python src/get_lightcurves.py crts --samples Graham2015,Charisi2016,Foustoul2025 --wait-hours 100 > data/lightcurves/crts_run.out 2>&1 &`
  As of 16:45 UTC: 32 Graham, 3 Charisi with CRTS data.
- ZTF Foustoul (239/257 done):
  `for k in 0 1 2; do nice .venv/bin/python src/get_lightcurves.py ztf --samples Foustoul2025 --shard $k/3 --wait-hours 24 & done`
  The other samples' ZTF is exploratory only (not needed for the pre-registration).
- simulate_null on 31 Charisi templates (K=200), 30/31 done: interim r0 = 0.06 (95% upper bound
  0.066), false survivors (lnB > 3) 0.55%.

## Results so far
- Charisi null calibration complete: 6200 fakes, r0 = 0.061 (95% CI 0.055-0.066), lnB>3 0.56%.
- Dress rehearsal (Charisi, mock DR4): null k=2/30, no survivors; signal H1 rejects, 7/19 pass H4.
- DHO sensitivity null sim done: 3100 fakes, r0 = 0.051 (95% upper 0.058), lnB>3 1.3%.
- PTF for Graham (fallback group forecast_alt) downloaded; provisional full freeze running
  -> results/frozen_provisional/ (not final; final freeze after CRTS passes, by ~28 Nov).
- Sensitivity modes + negative control in analyze_dr4.py (--sensitivity ...); figures.py.

## Provisional freeze (results/frozen_provisional, 18:00 UTC)
- forecast 63 scored (Graham-CRTS 32 + Charisi 31): 56 with mock power >= 0.5
- forecast_alt 43 (Graham via PTF): 39 with power >= 0.5;  gap 73 (Foustoul ZTF-only): 7
- 33 Graham with neither CRTS nor PTF (so far)
- Null sims for Graham templates running: 2 shards
  `nice .venv/bin/python src/simulate_null.py --k 200 --workers 1 --shard k/2 --frozen-dir results/frozen_provisional --groups forecast,forecast_alt`

## Design changes today (see notes/decisions.md)
- sigma_s fixed at 1 mag (broad); fake selection = power at the published period.
- Charisi mock power: 29/31 have power(p<0.05) >= 0.5.

## State at the stop point
- Freeze: results/frozen (67 forecast, 41 forecast_alt, 73 gap). Analysis code commit 9284e83.
- preregistration/ assembled by `src/make_prereg.py --code-hash 9284e83... --power-rule 0.5`,
  UNSTAGED. r0 = 0.0562, r0_used = 0.0596. Summary: notes/prereg_summary.md.
- Dress rehearsal on mock DR4 passed (results/dress_rehearsal_*).
- Gaia archive is down 25 Nov 12:00 CET -> 2 Dec; DR3 mean G cached; scripts tested on DR3.

## Next
1. On approval: `git add preregistration && git commit && git push` (before 2 Dec 12:00 CET).
   Recommended: merge research-plan into main first so it is on the default branch.
2. CRTS retries continue (detached). If more Graham objects arrive BEFORE approval: rerun
   `src/frozen.py`, drop changed templates from results/null_sim_drw*.csv, rerun
   `src/simulate_null.py --k 200 --workers 1 --groups forecast,forecast_alt --out results/null_sim_drw_refresh2.csv`,
   commit code, rerun make_prereg with the new hash.
3. After DR4 release and the public prereg commit: gaia_xmatch.py, gaia_epochs.py, analyze_dr4.py.

## Risks
- CRTS outage: Graham (111) is the largest forecast sample. If CRTS does not return before
  ~29 Nov, options are in the report (pre-register Charisi only plus a dated amendment for
  Graham before DR4; or use PTF/iPTF for Graham objects as a deviation).
