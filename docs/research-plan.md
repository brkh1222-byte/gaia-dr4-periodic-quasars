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

_(in progress)_

## A.3 Open questions and debates

_(in progress)_

## A.4 Is our question still unanswered?

_(in progress)_

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

## B.3 Analysis plan

_(in progress)_

## B.4 Positive, negative and null results

_(in progress)_

## B.5 Paper outline, figures and target journal

_(in progress)_

## B.6 Timeline and risks

_(in progress)_

## Suggested changes to PROPOSAL.md / TASKS.md

_(in progress)_
