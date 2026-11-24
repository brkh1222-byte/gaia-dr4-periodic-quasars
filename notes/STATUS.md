# STATUS — 2026-11-24 (evening)

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

## Running (background)
- Gaia DR3 mean G (gaia_source only) -> data/raw/dr3_mean_g.csv
- ZTF downloads (all 447) -> data/lightcurves/ztf/
- CRTS: server down since ~09:00 UTC 24 Nov; `get_lightcurves.py crts --wait-hours 96`
  probes every 15 min and downloads when it returns.

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
