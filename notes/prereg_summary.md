# Pre-registration in one page (needs your OK before 2 Dec 2026, 12:00 CET)

Full text: `preregistration/prereg.md` (assembled, not committed). Freeze of 28 Nov 2026;
analysis code commit 9284e83. Commands for the moment you say OK: `notes/prereg_commit_commands.md`.

## Three choices to confirm (made by Claude while you were away)

1. **Testability rule: at least 50% power.** A candidate counts in the headline test only if,
   were its sinusoid real, Gaia would give p < 0.05 at least half the time. Alternative: 80%
   (about a third fewer candidates, same error rate).
2. **Graham candidates without Catalina (CRTS) data are a separate, secondary group** that uses
   PTF data instead (41 objects). They are not in the headline test. Alternatives: drop them,
   or pool them into the headline test.
3. **Merge the work branch `research-plan` into `main` before the pre-registration commit**, so
   the timestamped commit sits on the repository's default branch.

## What is tested

Gaia DR4 photometry (Jul 2014 – Jan 2020) was never seen by the authors of the candidate
catalogues. Does each published period, with phase and amplitude worked out only from the
discovery data, predict that photometry better than red noise does?

**Frozen before DR4 (your option A).** No catalogue published a phase or reference epoch, and
Graham+2015 and Charisi+2016 gave no period errors. Each period is used exactly as printed.
Phase and amplitude come only from light-curve data taken before the paper's arXiv v1 date.
These values and both models' predicted light curves are in `frozen_params.csv` and
`frozen_predictions.csv`.

**Samples.** Headline: 67 candidates = 34 Graham+2015 (with CRTS data) + 33 Charisi+2016 (PTF),
all with Gaia G < 20.7 and at least 10 discovery nights; Gaia epochs after 27 Jul 2015 / 4 Apr
2016. Secondary: 41 Graham via PTF; 85 Foustoul+2025 in the gap before their ZTF data (low
power). Not scorable: 31 Graham objects with neither CRTS nor PTF data (the CRTS server has been
failing since 24 Nov), 5 with too few nights, and 158 Foustoul objects whose CRTS data could not
be retrieved; Liu+2019 and Chen+2024 (data or periods not public).
Nobody has looked at Gaia DR3 light curves of any candidate.

**Score.** ln B = how much better "sinusoid + red noise" predicts the Gaia points than "red
noise only", both conditioned on the discovery data. The Gaia zero point, amplitude ratio and
extra noise are free and treated identically in both models.

**Per-candidate p-value.** 1000 simulated Gaia data sets from the red-noise prediction at the
real DR4 times; p = fraction with ln B at least as large as observed.

**Headline test (H1).** Count candidates with p < 0.05 among the testable ones (60 of 67 on a
mock cadence). Excess if a binomial test gives P < 0.01 against a false-positive rate of 0.0596.
That rate was measured before release with 13,400 red-noise fake candidates run through the
whole pipeline (measured 0.056; 95% upper bound used).

**Individual candidates.** Survivor: ln B > 3 and Benjamini–Hochberg q < 0.05. "Ephemeris-specific"
if it also beats 95% of 200 phase-shifted forecasts. Never called "confirmed binary".

**Fraction.** 95% upper limit on the fraction f of candidates that are truly periodic.

**Checks reported alongside.** Shuffled periods; only epochs after the DR3 window; unbinned
Gaia data; period uncertainty included; wider amplitude prior; steeper red-noise spectrum.

**Practice run on fake Gaia data.** Red-noise data: 3 of 62 below p = 0.05 (3.1 expected), no
survivors, f < 0.05. All ephemerides true: excess detected, 38 survivors, f ≈ 0.94.
