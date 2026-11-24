# STATUS — 2026-11-25 00:45 local (2026-11-24 16:45 UTC)

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

## Design changes today (see notes/decisions.md)
- sigma_s fixed at 1 mag (broad); fake selection = power at the published period.
- Charisi mock power: 29/31 have power(p<0.05) >= 0.5.

## Next
1. `src/frozen.py` (needs DR3 + light curves) -> results/frozen/*; power on mock cadence.
2. `src/simulate_null.py --k 50 --workers 2` on forecast templates -> r0.
3. Dress rehearsal `src/analyze_dr4.py --mock null` and `--mock signal`.
4. Choose testability rule (open question), commit code, record hash,
   `src/make_prereg.py --code-hash <hash> --power-rule <x>`; DO NOT git add preregistration/.
5. notebooks/make_01_sample.py; finish paper Data/Methods numbers.

## Risks
- CRTS outage: Graham (111) is the largest forecast sample. If CRTS does not return before
  ~29 Nov, options are in the report (pre-register Charisi only plus a dated amendment for
  Graham before DR4; or use PTF/iPTF for Graham objects as a deviation).
