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
