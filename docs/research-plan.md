# Gaia DR4 Periodic Quasars — Research Understanding & Plan

Mirror of the Claude Doc https://claude.ai/artifact/XLmjm3H9pSoC948ju4T5xq (as of 2026-11-24). The doc is the working copy; this file is updated as sections are completed.

> Sections marked _(in progress)_ are not yet written.

## Summary

_(in progress)_

# Part A — What we're studying

## A.1 Core concepts

A periodic quasar is a quasar whose brightness seems to rise and fall on a regular clock of a few years; the hope is that the clock is two black holes orbiting each other. Our job is to check whether those clocks keep ticking when we look at years nobody used to find them.

**Quasars.** A quasar is a supermassive black hole (millions to billions of Suns) in a galaxy's centre, swallowing gas from a hot, glowing disk. The disk can outshine the whole galaxy, so we see a quasar as a point of light even billions of light-years away.

**Binary black holes (SMBHBs).** Most big galaxies have a central black hole, and galaxies merge. After a merger the two black holes sink to the centre and end up orbiting each other. When they are closer than about 0.01–0.1 parsec, the orbit takes only years. Kepler's third law gives the period:

```latex
P = 2\pi\sqrt{\frac{a^{3}}{G M}} \approx 3\ \text{yr}\times\left(\frac{a}{0.01\,\text{pc}}\right)^{3/2}\left(\frac{M}{10^{9}M_\odot}\right)^{-1/2}
```

We observe the period stretched by cosmic expansion, by a factor (1 + z), where z is the redshift. A quasar at z = 1 with a 3-year orbit shows a 6-year period on Earth.

**Why a binary could make light periodic.** Three mechanisms are proposed:

- **Modulated feeding.** Simulations show gas falling from a disk around both holes in lumps, once per orbit or at a few related periods. More gas means more light.
- **Doppler boosting.** The smaller black hole moves at a few percent of light speed. When it moves towards us its light is beamed and blue-shifted, so it looks brighter. This predicts a smooth sinusoid with a specific link between optical and ultraviolet amplitudes.
- **Self-lensing.** If the orbit is edge-on, one hole's gravity briefly magnifies the other, giving a short flare once per orbit.

**Red noise: the main villain.** Every quasar flickers randomly, even with one black hole. This flicker is "red": slow wanders are bigger than fast ones. The standard model is the **damped random walk (DRW)**. The brightness drifts randomly but is pulled back to its mean on a timescale τ, typically hundreds of days. Its correlation between two times falls off exponentially:

```latex
\operatorname{cov}\big(m(t_i), m(t_j)\big) = \sigma^{2}\exp\!\left(-\frac{|t_i - t_j|}{\tau}\right)
```

With only 1.5–3 cycles of data, a random wander often looks like a wave. Searching thousands of quasars for waves therefore finds many fake ones. This is the central reason the field is sceptical.

**The decisive test: predict the future (or an unseen past).** A real binary keeps its period. A red-noise fluke does not. So the cleanest test asks: does the published sine wave predict data that were never used to find it? Gaia DR4 covers mid-2014 to early 2020. That window falls between CRTS (to 2013) and ZTF (from 2018), so for CRTS- and PTF-era candidates it is fresh data.

**Gaia epoch photometry.** Gaia is an ESA spacecraft that spins and scans the whole sky. Each star or quasar is measured many times as the scan passes over it. Every measurement, with its time, is an "epoch". DR4 is the first release to publish these epochs for essentially all sources, not just pre-selected variables.

**Pre-registration.** We write down our exact predictions, models and decision rules, and timestamp them on GitHub before the data exist publicly. Nobody, including us, can then tune the analysis to get a preferred answer.

## A.2 Key literature

The story of the field since 2015 is "hundreds of candidates found, then almost all fail when more data arrive". Three threads connect the papers below: candidate searches, retests with longer baselines, and independent physics checks (gravitational waves, Doppler colours). Newest first. Entries marked † were confirmed from search results or other papers' reference lists, not from an opened abstract page. Check them before citing.

| Year | Paper | What it found | How it connects |
| --- | --- | --- | --- |
| 2026 | [Huijse et al., A&A 715, A92](https://arxiv.org/abs/2505.16884) | Blind search of 377,128 Gaia DR3 AGN. 13 pass p < 10⁻⁵, matching the false-positive rate; no reliable candidates. | Shows Gaia's ~3-year DR3 baseline is too short. Pipeline "can scale to DR4": the likeliest fast DR4 search. |
| 2026 | [El-Badry, Hogg & Rix, PASP 138, 024102](https://arxiv.org/abs/2509.10601) | Gaia DR3 sinusoids fail to predict ZTF for all 116 testable objects. Stable sinusoids: ≤ a few per 10⁶ AGN. | Closest in spirit to our test (out-of-sample prediction), but Gaia → ZTF and Gaia-selected only. |
| 2026 | [Lin, Charisi & Haiman, ApJ 997, 316](https://arxiv.org/abs/2505.14778) | Lomb–Scargle misses most non-sinusoidal (sawtooth) binary signals. | Warns that a sinusoid-only test can miss real binaries. Limits what a negative result means. |
| 2026 | [Agarwal et al. (NANOGrav), ApJL 998, L11](https://arxiv.org/abs/2508.16534) | Targeted gravitational-wave search of 114 candidate AGN: no individual signal. | Independent physics check; no detection. |
| 2025 | [Foustoul et al., A&A 699, A55](https://arxiv.org/abs/2505.06656) | CRTS + ZTF (~15 yr): 36 strong and 221 weaker candidates. Of 89 Graham candidates with ZTF, 58 look like red noise and 26 consistent. | Closest check of CRTS ephemerides, but by refitting, not fixed prediction. Also a source of new candidates. |
| 2025 | [Casey-Clyde et al., ApJ](https://arxiv.org/abs/2405.19406) † | The NANOGrav background limits how many CRTS quasars can host binaries. | Physics prior on how many survivors to expect. |
| 2024 | [Robnik et al., MNRAS](https://arxiv.org/abs/2407.17565) † | Re-analysis of Charisi's 33 PTF candidates with a better null: none significant. | The PTF sample is already weakened. |
| 2024 | [Y.-J. Chen et al., MNRAS 527, 12154](https://arxiv.org/abs/2206.11497) | ZTF: 184 → 86 candidates; 3 keep their periodicity with archival data. | ZTF-era candidates; for them Gaia is a retrodiction. |
| 2024 | [Davis et al., ApJ](https://arxiv.org/abs/2311.10851) † | Periodicity-only selection yields samples dominated by false positives (~60% among luminous quasars). | Sets a low prior. |
| 2023 | [NANOGrav 15-yr, ApJL 951, L8](https://arxiv.org/abs/2306.16213) | Evidence for a nanohertz gravitational-wave background. | Main motivation: the background is thought to come from SMBHBs. |
| 2023 | [D'Orazio & Charisi, review](https://arxiv.org/abs/2310.16896) | Standard review of binary signatures. No close binary confirmed. | Best single background read for the student. |
| 2022 | [Witt et al., ApJ](https://arxiv.org/abs/2110.07465) † | Bayesian periodic + DRW detection for CRTS- and LSST-like data. | Template for our model-comparison approach. |
| 2022 | [Casey-Clyde et al., ApJ 924, 93](https://arxiv.org/abs/2107.11390) | Binary population model tied to quasars. | Expected numbers. |
| 2020 | [Zhu & Thrane, ApJ 900, 117](https://arxiv.org/abs/2004.10944) | PG 1302-102: quasi-periodic model slightly beats a sinusoid. | Real binaries may not be strict sinusoids. |
| 2019 | [Liu et al., ApJ 884, 36](https://arxiv.org/abs/1906.08315) | PS1 MDS: 26 → 1 after extended monitoring. | Shows how candidates decay with longer baselines. |
| 2018 | [Sesana, Haiman, Kocsis & Kelley, ApJ 856, 42](https://arxiv.org/abs/1703.10611) | If CRTS/PTF candidates were real, the gravitational-wave background would exceed PTA limits. False positives favoured (2.3σ CRTS, 3.6σ PTF). | Independent evidence most candidates are false. |
| 2018 | [Liu, Gezari & Miller, ApJL](https://arxiv.org/abs/1803.05448) | ASAS-SN extension weakens PG 1302-102's periodicity. | Early extended-baseline test. |
| 2018 | [Charisi et al., MNRAS 476, 4617](https://arxiv.org/abs/1801.06189) † | At most about a third of candidates fit the Doppler-boost colour test. | Independent physics check. |
| 2016 | [Vaughan et al., MNRAS 461, 3145](https://arxiv.org/abs/1606.02620) | Red noise readily fakes 2–3 cycles; PG 1302-102 fits red noise. | Origin of the "≥ 5 cycles" rule of thumb. |
| 2016 | [Charisi et al., MNRAS 463, 2145](https://arxiv.org/abs/1604.01020) | PTF: 33 candidates. | Input catalogue. |
| 2016 | [Liu et al., ApJ 833, 6](https://arxiv.org/abs/1609.09503) | PS1 MD09: 3 candidates, 0 survive. | Early extended-baseline failure. |
| 2015 | [Graham et al., MNRAS 453, 1562](https://arxiv.org/abs/1507.07603) | CRTS: 111 candidates from ~243,500 quasars. | Main input catalogue. |
| 2015 | [Graham et al., Nature](https://arxiv.org/abs/1501.01375) † | PG 1302-102, period ~1884 d. | The famous single candidate. |
| 2015 | [D'Orazio, Haiman & Schiminovich, Nature](https://arxiv.org/abs/1509.04301) † | Doppler-boost model for PG 1302-102. | Physical basis for a sinusoid. |
| 2009 | Kelly, Bechtold & Siemiginowska, ApJ 698, 895 | Quasar variability is well described by a damped random walk. | Our null model. |

## A.3 Open questions and debates

The field agrees that most candidates are red noise. It disagrees on what a real binary's light curve looks like, and therefore on how strict a test should be.

1. **Which red-noise model is the right null?** DRW, power-law, CARMA and newer "learned" null models give different significances ([Vaughan 2016](https://arxiv.org/abs/1606.02620); [Zhu & Thrane 2020](https://arxiv.org/abs/2004.10944); [Sellam et al. 2026](https://arxiv.org/abs/2611.15040)). *For us:* fix the null model and its calibration in advance, and run a sensitivity check with a steeper-spectrum null.
2. **How many cycles are enough?** The rule of thumb is ≥ 5 ([Vaughan 2016](https://arxiv.org/abs/1606.02620)). DR4's 5.5 years covers only ~1–4 cycles for 500–2000-day periods. *For us:* test whether the predicted phase and amplitude hold, rather than re-detecting a period in Gaia alone.
3. **Strictly periodic or quasi-periodic?** Real binaries may drift in phase, change amplitude, or have sawtooth shapes ([Zhu & Thrane 2020](https://arxiv.org/abs/2004.10944); [Lin et al. 2026](https://arxiv.org/abs/2505.14778)). *For us:* the allowed phase tolerance must be fixed before the data, or the test cannot fail.
4. **Do the candidates fit the gravitational-wave background?** The candidate lists imply more binaries than pulsar timing allows ([Sesana et al. 2018](https://arxiv.org/abs/1703.10611); [Casey-Clyde et al. 2025](https://arxiv.org/abs/2405.19406)). The NANOGrav targeted search found nothing ([Agarwal et al. 2026](https://arxiv.org/abs/2508.16534)). *For us:* our upper limit on f feeds straight into this debate.
5. **Why are stable sinusoids so rare?** Gaia data allow ≤ a few per 10⁶ AGN ([El-Badry et al. 2026](https://arxiv.org/abs/2509.10601)), while population models predict dozens of binaries in Gaia. Either binary light curves are not sinusoidal or binaries are very rare.
6. **Do Doppler and spectroscopic tests agree?** At most about a third of candidates pass the UV/optical colour test ([Charisi et al. 2018](https://arxiv.org/abs/1801.06189)). Spectroscopic tests of PG 1302-102 are inconclusive.
7. **Can predictions be trusted after the fact?** Predicted events have failed or been disputed: the OJ 287 outburst expected in Oct 2022 was not seen ([Komossa et al. 2023](https://arxiv.org/abs/2302.11646)). The predicted merger of SDSS J1430+2303 did not occur ([Jiang et al. 2022](https://arxiv.org/abs/2201.11633); [Dotti et al. 2023](https://arxiv.org/abs/2205.06275)). No sample-level pre-registration exists. *For us:* this is the methodological gap we fill.
8. **Gaia-specific systematics.** Gaia's cadence is sparse and irregular. Scan direction and host-galaxy light affect G-band photometry of extended sources ([Huijse et al. 2026](https://arxiv.org/abs/2505.16884)). *For us:* pre-specify the cross-band calibration and quality cuts.

## A.4 Is our question still unanswered?

**Yes, as of 24 Nov 2026, but the gap is narrower than PROPOSAL.md says, and two framing problems must be fixed before the pre-registration commit.** We found no paper, and no announced DR4 paper, that tests fixed published ground-survey ephemerides against Gaia photometry, or that pre-registers predictions for a sample of candidates. The search covered arXiv, ADS-indexed journals and ESA's [DR4 papers list](https://www.cosmos.esa.int/web/gaia/dr4-papers). It was web-based and not exhaustive.

**Closest prior work and how we differ.**

| Paper | What it did | How our test differs |
| --- | --- | --- |
| [El-Badry, Hogg & Rix 2026](https://arxiv.org/abs/2509.10601) | Asked whether Gaia DR3 sinusoids predict later ZTF data. All 116 failed. | Opposite direction (ground ephemeris → Gaia). Tests published survey catalogues. Pre-registered. |
| [Foustoul et al. 2025](https://arxiv.org/abs/2505.06656) | Compared 89 Graham candidates with ZTF: 58 red noise, 26 consistent, 5 unclear. | Fixed published ephemerides, not refits. Independent instrument. Fills the 2014–2018 gap. Pre-registered. |
| [Robnik et al. 2024](https://arxiv.org/abs/2407.17565) | Re-analysed Charisi's 33 PTF candidates: none significant. | New data rather than re-analysis of old data. |
| [Huijse et al. 2026](https://arxiv.org/abs/2505.16884) | Blind Gaia DR3 search for new periodic AGN. | We test published claims; they search for new ones. They are the group most likely to publish a DR4 search quickly. |

**Framing problem 1: part of the "blind" window is already public.** Gaia DR3 (2022) published epoch photometry for about a million variable-quasar candidates. It covers roughly Jul 2014 to May 2017, about the first 34 months of the DR4 window. Huijse et al. used it. Some of our candidates may already have public Gaia light curves.

- The CRTS and PTF ephemerides (2015–2016) were published before any Gaia quasar light curves existed, so the test of *their* predictions stays out-of-sample.
- What DR3 threatens is *our* blindness as analysts. The pre-registration must state that nobody on the project has inspected DR3 epoch photometry for any candidate, and must forbid doing so.
- Add a robustness test that uses only the DR4-new interval (~May 2017 – Jan 2020).
- The DR3 end date is from secondary sources. Confirm it against the ESA DR3 documentation before committing.

**Framing problem 2: for some catalogues DR4 is not "the future".** Chen 2024 and Foustoul 2025 fitted ZTF data starting in 2018, inside the DR4 window. Liu 2019's follow-up likely ran into 2015–2018. For these the honest phrase is "independent-instrument test of fixed published ephemerides". Each catalogue's group label (forecast or retrodiction) and its usable Gaia epochs must be pre-registered.

**Scoop risk.** ESA's DR4 paper list includes a second variable-AGN catalogue and lensed-quasar time delays, but no periodic-AGN or SMBHB paper. The main external risk is a fast DR4 periodicity search by the Huijse group or the El-Badry/Hogg/Rix group. Neither would test published ephemerides. A dated pre-registration makes priority on the *method* unambiguous even if results land in the same month.

## A.5 Standard methods, datasets, tools and pitfalls

_(in progress)_

## A.6 Verified references (citation check of PROPOSAL.md)

Two references in PROPOSAL.md are seriously wrong and five more need fixes or missing details. Each was checked on its arXiv abstract page and against Crossref metadata on 24 Nov 2026.

| # | Reference as verified | Verdict | What to fix in PROPOSAL.md |
| --- | --- | --- | --- |
| 1 | Graham, Djorgovski, Stern et al. 2015, MNRAS 453, 1562. [arXiv:1507.07603](https://arxiv.org/abs/1507.07603), doi:10.1093/mnras/stv1726 | Verified | None. ~243,500 quasars searched, 111 candidates. |
| 2 | Charisi, Bartos, Haiman et al. 2016, MNRAS 463, 2145. [arXiv:1604.01020](https://arxiv.org/abs/1604.01020), doi:10.1093/mnras/stw1838 | Verified | None. 35,383 quasars; 50 significant vs DRW, 33 still significant with iPTF + CRTS. |
| 3a | Liu, Gezari, Ayers et al. 2019, ApJ 884, 36. [arXiv:1906.08315](https://arxiv.org/abs/1906.08315), doi:10.3847/1538-4357/ab40cb | Verified | Note only 1 of 26 candidates stays significant vs DRW after extended monitoring. |
| 3b | Liu, Gezari, Burgett et al. 2016, ApJ 833, 6. [arXiv:1609.09503](https://arxiv.org/abs/1609.09503), doi:10.3847/0004-637X/833/1/6 | Verified | Add journal. 3 candidates; none survived an extended baseline. |
| 4 | **Y.-J.** Chen, Zhai, Liu et al. **2024**, MNRAS 527, 12154. [arXiv:2206.11497](https://arxiv.org/abs/2206.11497), doi:10.1093/mnras/stad3981 | Minor error | First author is Y.-J. Chen (not Y.-C.). One paper, not "2022/24". "3 final" means 3 of 86 kept their periodicity with archival extension. |
| 5 | Foustoul, Webb, Mignon-Risse, Kammoun, Volonteri & Dong-Páez 2025, A&A 699, A55. [arXiv:2505.06656](https://arxiv.org/abs/2505.06656), doi:10.1051/0004-6361/202452818 | Minor error | Add authors and journal. Numbers are right: 36 strong + 221 weaker candidates; 58 earlier CRTS candidates attributed to red noise. |
| 6 | Huijse, Davelaar, De Ridder, Jannsen & Aerts 2026, A&A 715, A92. [arXiv:2505.16884](https://arxiv.org/abs/2505.16884), doi:10.1051/0004-6361/202660885 | **Major error** | v2 finds 13 sources at p < 10⁻⁵, consistent with false positives: "no reliable periodic SMBHB candidates". Do not write "181 → 13 candidates". The paper says only that its pipeline can scale to DR4, not that a blind DR4 search is planned. |
| 7 | El-Badry, Hogg & Rix 2026, PASP 138, 024102. [arXiv:2509.10601](https://arxiv.org/abs/2509.10601), doi:10.1088/1538-3873/ae3f56 | Verified | Add volume and page. Tested Huijse v1's 181 candidates. All 116 with ZTF failed. Strictly periodic AGN are at most a few per 10⁶. |
| 8 | Lin, Charisi & Haiman **2026**, ApJ 997, 316. [arXiv:2505.14778](https://arxiv.org/abs/2505.14778), doi:10.3847/1538-4357/ae29a7 | Minor error | Year is 2026. Lomb–Scargle recovers only ~45% (PTF-like) of injected sinusoids and far fewer sawtooth signals. |
| 9 | Vaughan, Uttley, Markowitz et al. 2016, MNRAS 461, 3145. [arXiv:1606.02620](https://arxiv.org/abs/1606.02620), doi:10.1093/mnras/stw1412 | Verified | Add arXiv ID. It is about PG 1302-102 and red-noise false periods. |
| 10 | Kelly, Bechtold & Siemiginowska 2009, ApJ 698, 895. doi:10.1088/0004-637X/698/1/895 | Verified | None. |
| 11 | Zhu & Thrane 2020, ApJ 900, 117. [arXiv:2004.10944](https://arxiv.org/abs/2004.10944), doi:10.3847/1538-4357/abac5a | **Major error** | The proposal calls it a red-noise test of the Graham 2015 sample. It is a Bayesian study of **PG 1302-102 only**. |
| 12 | Gaia DR4: release 2 Dec 2026; data 25 Jul 2014 – 20 Jan 2020 (66 months). [ESA DR4 page](https://www.cosmos.esa.int/web/gaia/data-release-4), [DR4 content](https://www.cosmos.esa.int/web/gaia/dr4) | Minor error | ESA gives no "~2 billion with epoch photometry" count. Epoch data cover all ~2.8 billion processed sources (`epoch_photometry`: ~1.2 × 10¹¹ rows, 7.2 TB). |
| 13 | NANOGrav: Agazie et al. 2023, ApJL 951, L8. [arXiv:2306.16213](https://arxiv.org/abs/2306.16213), doi:10.3847/2041-8213/acdac6 | Verified | Add as the PTA citation in §2. |

**DR4 facts for the plan.** DR4 includes a `vari_agn` table (~1.97 million rows) and a `qso_candidates` table (~5.4 million rows). None of the ~87 planned DPAC DR4 papers ([list](https://www.cosmos.esa.int/web/gaia/dr4-papers)) covers periodic AGN or SMBHBs. Astroquery's development code lists DR4 retrieval types, up to 5,000 IDs per call. The exact string `data_release='Gaia DR4'` must be confirmed on release day ([astroquery source](https://raw.githubusercontent.com/astropy/astroquery/main/astroquery/gaia/core.py)).

# Part B — How we'll study it

## B.1 Refined research question and hypotheses

**Research question.** For periodic-quasar candidates whose discovery data end before mid-2014, does the published sinusoidal ephemeris predict their Gaia DR4 G-band light curve (Jul 2014 – Jan 2020) better than a damped-random-walk model conditioned on the same discovery data? And what fraction of candidates can be true periodic sources?

Three changes sharpen the proposal's question:

1. **Per-candidate training cut-off.** Each candidate is scored only on Gaia epochs after the last date in its own discovery data. This keeps every test out-of-sample. ZTF-era candidates become a separate "retrodiction" group.
2. **A measured quantity, not just a verdict.** We estimate the true-binary fraction f (with an upper limit if no excess). A null result then still yields a number.
3. **A power criterion before scoring.** Simulations decide, before release, which candidates Gaia can actually test. Untestable ones are reported but excluded from the primary test.

**Hypotheses (pre-registered).**

| ID | Statement | Test |
| --- | --- | --- |
| H0 (null) | Published ephemerides carry no predictive information beyond red noise. The log Bayes factors (lnB) of testable candidates follow the DRW-null distribution. | lnB distribution vs simulated null |
| H1 (population) | More candidates than expected under the null have lnB above the pre-set threshold. | One-sided binomial test of the count |
| H2 (individual) | Specific candidates show lnB above threshold at a false-discovery rate below 5%. | Null-calibrated p-values, Benjamini–Hochberg |
| H3 (fraction) | The fraction f of testable candidates that are true periodic sources is constrained. | Mixture of null and injected-signal lnB distributions |
| H4 (ephemeris specificity, control) | Any signal depends on the published phase and period, not on any sinusoid. | Same score with scrambled phases and perturbed periods |

**Expected outcome.** Most candidates fail, given recent ZTF-based retests. The informative quantities are the upper limit on f and any individual survivors.

## B.2 Study design and methods

_(in progress)_

## B.3 Analysis plan: statistics, controls, decision rules

Every candidate gets one number, the predictive log Bayes factor lnB. Simulations, not textbook thresholds, decide what that number means. All rules below are frozen in the pre-registration.

**1. Per-candidate score.** Let D_old be the discovery data and D_new the Gaia epochs after that candidate's cut-off. Each model is conditioned on D_old and asked how probable D_new was:

```latex
\ln B = \ln p(D_\text{new} \mid D_\text{old}, \mathcal{M}_\text{sin+DRW}) - \ln p(D_\text{new} \mid D_\text{old}, \mathcal{M}_\text{DRW})
```

Both predictive densities are Gaussian-process conditionals computed with celerite2. lnB > 0 means the sinusoid predicted Gaia better. Nuisance parameters are fixed in advance and marginalised the same way in both models:

- a Gaia–survey magnitude offset (analytic, broad Gaussian prior);
- a band amplitude ratio G/V (grid over a pre-set log-normal prior);
- a Gaia white-noise jitter (grid, prior centred on Gaia's documented per-transit precision at the candidate's G magnitude, fixed before release).

**2. Null calibration (what lnB looks like if the candidate is red noise).** For each candidate, simulate 1,000 DRW light curves using its own discovery cadence, errors and fitted DRW parameters. Keep only the simulations that pass a simplified copy of the discovery selection: periodogram peak in the searched range, at least 1.5 cycles, and amplitude at least as large as published. Fit their "ephemerides", extend each light curve onto Gaia's cadence and noise, and score them. This gives a null distribution and a p-value pᵢ for every candidate.

**3. Injection (what lnB looks like if the candidate is real).** Simulate the published sinusoid plus DRW on the same cadences and score it. This gives the signal distribution and each candidate's power.

**4. Testability cut (decided before release).** A candidate is *testable* if injected signals exceed the threshold (lnB > 3) in at least 80% of simulations. Untestable candidates are listed and excluded from H1–H3.

**5. Pre-registered tests.**

| Hypothesis | Statistic | Decision rule |
| --- | --- | --- |
| H1 population excess | k = number of testable candidates with pᵢ < 0.05 | Reject H0 if the one-sided binomial p (n, 0.05) < 0.01 |
| H2 individual survivors | pᵢ for each testable candidate | Survivor if lnB > 3 and Benjamini–Hochberg q < 0.05 |
| H3 binary fraction f | Likelihood ∏ᵢ [f·p_inj,ᵢ(lnBᵢ) + (1−f)·p_null,ᵢ(lnBᵢ)], flat prior on f | Report posterior median and 95% upper limit |
| H4 ephemeris specificity | lnB with published phase vs 200 random phases, and period ±10% | A survivor counts only if its published-phase lnB beats 95% of scrambled versions |

**6. Controls.**

- **Phase scrambling and period perturbation** (H4) show whether the specific published ephemeris matters, rather than any wave.
- **Negative control:** swap ephemerides between candidates at random. The H1 count should fall to the null rate.
- **Positive control:** inject fake binaries into the real Gaia cadence of non-candidate quasars. The pipeline must recover them at the rate step 3 predicts.
- **Dress rehearsal:** run the whole pipeline end to end on mock DR4 files before release. Freeze the code by commit hash.

**7. Pre-registered sensitivity analyses (reported, not used for the headline).**

- A steeper-spectrum null (CARMA(2,1)/damped harmonic oscillator) instead of the DRW.
- Only the DR4-new interval (after DR3's May 2017 end).
- Unbinned versus visit-binned Gaia data.
- A wider amplitude-ratio prior.
- A quasi-periodic variant allowing slow phase drift.

**8. Exploratory (labelled as such in the paper).** Joint CRTS + Gaia + ZTF fits for any survivors. Comparison with Foustoul's ZTF verdicts for the same objects.

## B.4 Positive, negative and null results

All three outcomes are publishable, because the pre-registration and the power analysis are fixed before we see the data. A reader can trust a "no" as much as a "yes".

| Outcome | What we would see | What it means | Why it is publishable |
| --- | --- | --- | --- |
| **Positive** | H1 rejects the null (binomial p < 0.01), or at least one H2 survivor passes the H4 specificity check | Some published ephemerides correctly predicted 2014–2020 Gaia data they never saw | First pre-registered, out-of-sample confirmation of SMBHB ephemerides. Survivors become top targets for spectroscopy, pulsar-timing searches and Rubin/LSST monitoring. A strong RNAAS note, likely followed by an ApJ/MNRAS paper with the mentor. |
| **Negative** | Testable candidates' lnB values match the null, and the H3 upper limit on f is well below 1 | The ephemerides have no predictive power beyond red noise | Converts a qualitative "probably red noise" into a number (f < X at 95%) from independent data. Adds an independent instrument and the 2014–2018 gap to Foustoul's and El-Badry's results. Directly relevant to the debate over whether candidates fit the PTA background. |
| **Null / inconclusive** | Too few candidates are testable, or the H3 posterior on f stays broad | Gaia's cadence and depth cannot decide for these candidates | Still a useful methods note: a public, pre-registered pipeline and a power map of which candidate types Gaia can test. The frozen predictions can be re-scored against Rubin/LSST or later Gaia releases. Weaker; RNAAS-sized. |

**Guard against spin.** The headline is the pre-registered population result (H1 and H3). Individual survivors are reported with their q-values and H4 outcome, never as "confirmed binaries". Sensitivity and exploratory results are labelled as such. Any deviation from the plan goes into a dated amendment that explains why.

## B.5 Paper outline, figures and target journal

_(in progress)_

## B.6 Timeline and risks

_(in progress)_

## Suggested changes to PROPOSAL.md / TASKS.md

_(in progress)_
