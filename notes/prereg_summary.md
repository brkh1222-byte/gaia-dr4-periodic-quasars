# Pre-registration in one page (for the researcher's approval)

Status: DRAFT awaiting approval. Numbers are from the freeze of 24 Nov 2026 (analysis code commit
9284e83). If later CRTS retries add Graham objects before approval, the freeze is rerun and these
numbers updated. Full text: `preregistration/prereg.md`
(not yet committed). Decisions taken while you were away: `notes/decisions.md`.

**Question.** Gaia DR4 photometry (Jul 2014 – Jan 2020) was never seen by the authors of the
candidate catalogues. Does each catalogue's published period, with phase and amplitude worked out
only from its discovery data, predict that photometry better than red noise does?

**What is frozen (your option A).** No catalogue published a phase or reference epoch, and
Graham+2015 and Charisi+2016 published no period errors. So each period is used exactly as
printed. Phase and amplitude come from the light-curve data taken before the paper's arXiv v1
date (a Gaussian-process fit with a damped-random-walk noise model and a broad prior on the
sinusoid). These values, and both models' predicted light curves every 10 days, are frozen in
`frozen_params.csv` / `frozen_predictions.csv` before DR4 exists.

**Samples and cuts.**
- Primary: Graham+2015 CRTS (111) and Charisi+2016 PTF (final 33). Gaia epochs after each paper's
  arXiv v1 date (27 Jul 2015, 4 Apr 2016). DR3 mean G < 20.7; ≥ 10 discovery nights.
  67 scored: 34 Graham with CRTS data + 33 Charisi.
- Secondary, reported separately: Graham objects whose CRTS data could not be retrieved
  (the CRTS server is intermittent) are conditioned on PTF data instead (41). Foustoul+2025
  objects are tested in the gap before their ZTF data (73 scored so far, low power: 7 testable).
  36 Graham objects have neither CRTS nor PTF data yet and cannot be scored.
- Listed but not scored: Liu+2019 (light curves not public), Chen+2024 (periods not tabulated).
- Nobody has looked at Gaia DR3 epoch photometry for any candidate. A robustness test uses only
  DR4 epochs after the DR3 window (after 28 May 2017).

**Scoring.** For each candidate: ln B = log-probability of the Gaia data under "sinusoid +
red noise" minus that under "red noise only", both conditioned on the discovery data. The Gaia
zero point, the Gaia/survey amplitude ratio and extra Gaia noise are free nuisance parameters,
treated identically in both models. ln B > 0 means the published sinusoid predicted Gaia better.

**Per-candidate p-value.** Simulate 1000 Gaia data sets from the red-noise prediction (at the real
DR4 times) and count how often ln B is at least the observed value.

**Testable.** At least 50% power, i.e. if the sinusoid were real, p < 0.05 in at least half of the
simulations, computed from DR4 times and errors only. 60 of 67 primary candidates on a mock cadence (37 of 41 in the PTF-fallback group).

**Population test (H1, the headline).** k = number of testable primary candidates with p < 0.05.
"Excess" if a binomial test against rate r0 gives P < 0.01. r0 is measured before release
by simulating red-noise "fake candidates" that look as periodic as the real ones and running
the whole pipeline: r0 = 0.056 (13,400 fakes on the 67 primary templates; 95% upper bound 0.0596, which is
the rate used). For 10 strongly periodic Graham templates no fake matched in 200,000 trials, so
the 200 most periodic fakes were used (flagged "relaxed").

**Individual candidates (H2/H4).** Survivor: ln B > 3 and Benjamini–Hochberg q < 0.05.
"Ephemeris-specific" if it also beats 95% of 200 phase-shifted forecasts. Never called "confirmed".

**Fraction (H3).** Posterior on the fraction f of true periodic sources; report the 95% upper limit.

**Controls and sensitivity (not headline).** Shuffled periods (should give ~n·r0); DR4-only
epochs; unbinned Gaia transits; period marginalised over its error (5% if none); wider amplitude
ratio prior; steeper-spectrum (DHO) red noise (r0 = 0.051 for PTF).

**Dress rehearsal (mock DR4, 67 primary candidates, 62 testable).** With red-noise data: k = 3
(3.1 expected), no survivors, f < 0.05. With every ephemeris true: H1 rejected, 38 survivors (16
ephemeris-specific), f ≈ 0.94. Shuffled-period control on red-noise data: k = 2 of 45.
