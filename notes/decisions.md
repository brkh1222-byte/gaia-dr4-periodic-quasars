# Decisions log

## 2026-11-24 — What "published ephemeris" means (researcher: option A)
No source paper (Graham+2015, Charisi+2016, Liu+2019, Chen+2024, Foustoul+2025) publishes a
phase or reference epoch (T0). Graham+2015 and Charisi+2016 publish no period uncertainty.
Decision: freeze each published period exactly as printed; derive phase and amplitude only from
discovery data obtained before each paper's arXiv v1 date (GP posterior of the sinusoid
coefficients); freeze them in preregistration/frozen_params.csv before any DR4 data exist.
Say so explicitly in the pre-registration and in the paper.
Period marginalisation over the published error (or 5% when none) is a sensitivity analysis only.

## 2026-11-24 — Calibration design (Claude, for approval with the pre-registration)
- Per-candidate p-value: conditional null, Gaia data simulated from model D's predictive given
  the real discovery data (selection of the candidate is then ignorable).
- The H1 binomial rate r0 comes from the end-to-end simulation of selected DRW "fake candidates"
  (simulate_null.py), computed before release: r0_used = max(0.05, 95% upper bound).
- Positions: catalogue-verified (DR16Q/Milliquas) positions are used for all downloads and
  cross-matches; published positions are kept in the catalogue.

## 2026-11-24 — Sinusoid prior width (Claude; for approval with the pre-registration)
Fitting sigma_s (prior sd of the sinusoid coefficients) by evidence maximisation drove it to ~0
for about half of the Charisi candidates. The Occam penalty with 1.5-3 cycles let the DRW absorb
the sinusoid, i.e. the fit re-decided on the discovery data whether the sinusoid exists. Under
option A the published sinusoid is taken as given, so sigma_s is fixed at a broad 1 mag and the
phase and amplitude posterior is set by the discovery data. Effect on Charisi (mock cadence):
power(p<0.05) >= 0.5 for 29 of 31 instead of ~16 of 31.

## 2026-11-24 — Fake-candidate selection in the null simulation
A period-window peak search was rejected: with red noise the evidence periodogram rises toward
long periods, so the maximum sat at the window edge for nearly every fake and for the real
candidates too. Fakes are now selected by their red-noise-weighted power at the published
period being at least the real candidate's.
