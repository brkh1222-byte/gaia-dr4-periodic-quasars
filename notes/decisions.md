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

## 2026-11-24 — Open choices decided by Claude (researcher away; listed in the report for review)
1. **Testability rule: power(p < 0.05) >= 0.5**, computed from the real DR4 timestamps and
   errors only (1000 model-S draws). Reasons: the H1 statistic is the count of p < 0.05, so this
   is power for the statistic actually used; the binomial test keeps its size whatever rule is
   chosen (the rule only affects power); >= 0.8 would drop about a third of informative
   candidates (forecast mock: 56 vs 35 of 63); the research plan's P(lnB > 3) >= 0.8 keeps only
   20 of 63.
2. **Graham objects without CRTS data at the freeze:** use PTF/iPTF R data before the cutoff as
   the conditioning data, in a separate group `forecast_alt`, reported separately and NOT part of
   the primary H1–H3 test (the primary test uses only discovery-survey data: Graham-CRTS and
   Charisi-PTF). Objects with neither are listed as untestable. Reasons: stays faithful to option
   A for the headline while not discarding blind information. A pooled forecast + forecast_alt
   H1 is reported as a secondary analysis with the same r0_used.
3. **Survivor tiers:** "survivor" = ln B > 3 and BH q < 0.05 (H2). "Ephemeris-specific survivor"
   = survivor that also passes H4 (beats 95% of 200 phase-shifted forecasts). Reason: in the dress
   rehearsal with every ephemeris true, only 7 of 19 candidates with ln B > 3 passed H4. H4 tests
   phase specificity, which a real binary with a poorly constrained discovery phase can fail, so
   making it a gate would mostly lose true signals. Both tiers are reported; neither is called
   a "confirmed binary".
4. **Commit trailer:** commits up to 2026-11-24 ~18:30 UTC carry "Claude Opus 5.5", which was the
   attribution configured in the harness at the time. From the resume at ~19:00 UTC the harness
   attribution is "Claude Fable 5.2", and commits use that trailer.

## 2026-11-24 — Gaia archive downtime
ESA: the Gaia archive is down from 25 Nov 12:00 CET until the DR4 release (2 Dec 12:00 CET).
DR3 mean G for all candidates was cached on 24 Nov (data/raw/dr3_mean_g.csv). Release-day
scripts were tested against DR3 on 24 Nov: cross-match 179/179; epoch download + cleaning on
non-candidate control stars (data_structure must be DATAMODEL_STANDARD; INDIVIDUAL is deprecated).
