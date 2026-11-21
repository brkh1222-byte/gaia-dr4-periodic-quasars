# Research Proposal

## A Blind Test of Periodic-Quasar Candidates with Gaia DR4: Do Published Binary-Black-Hole Ephemerides Predict 2014–2020 Light Curves Better Than Red Noise?

**Author:** Berkh (independent student researcher)
**Mentor / prospective co-author:** [mentor name] (cosmology & black holes)
**Field:** Time-domain extragalactic astrophysics / supermassive black hole binaries
**Planned duration:** about 3 weeks, from Gaia DR4 release (**2 Dec 2026**) to a submitted RNAAS note, plus about 1 week of preparation before release
**Status:** Proposal, Nov 2026

---

## 1. Summary

Surveys such as CRTS, PTF and ZTF have reported hundreds of quasars whose brightness seems to oscillate with periods of a few years. These are proposed as **supermassive black hole binaries (SMBHBs)** at sub-parsec separation. Most candidates show only 1.5–3 cycles, and red-noise quasar variability readily mimics such periods. The field is now sceptical:
- El-Badry, Hogg & Rix (2026) showed that all 116 Gaia DR3 candidates with ZTF data fail.
- Foustoul et al. (2025) judged 58 earlier candidates to be consistent with red noise.

**Gaia DR4 (released 2 Dec 2026)** publishes G/BP/RP epoch photometry for about 2 billion sources, covering **25 Jul 2014 – 20 Jan 2020**. That window sits **in the gap between CRTS (ended 2013) and ZTF (started 2018)**. It is out-of-sample for every CRTS- and PTF-era candidate.

This project performs a **pre-registered blind test**:
- **Before** DR4 is released, the published sinusoidal ephemerides of each candidate are frozen and their predicted 2014–2020 light curves committed to this public repository.
- **After** release, each prediction is compared with the DR4 data. The test asks whether the sinusoid predicts better than a damped-random-walk (red-noise) model fitted to the same pre-Gaia data.

## 2. Background

- **SMBHBs** are expected from galaxy mergers and are the likely source of the nanohertz gravitational-wave background detected by pulsar timing arrays (2023). Confirming individual binaries is a major open problem.
- **Periodicity searches** mostly rely on short baselines. The damped random walk (DRW) describes quasar variability well and produces false periods at a known rate (Vaughan et al. 2016).
- **Out-of-sample prediction** is the cleanest test. Pre-registration before the data exist removes any possibility of tuning to the answer. With DR4 releasing on a known date, this is a rare opportunity for a truly blind test.

## 3. Prior work and the gap

| Sample | Size | arXiv | Already tested by |
|---|---|---|---|
| Graham et al. 2015, CRTS | 111 | 1507.07603 | Red-noise critiques (2004.10944); ZTF (Foustoul 2025) |
| Charisi et al. 2016, PTF | 33 | 1604.01020 | Partly Foustoul 2025 |
| Liu et al. 2016/2019, PS1 MDS | few | 1609.09503, 1906.08315 | Mostly too faint for Gaia |
| Chen et al. 2022/24, ZTF | 3 final | 2206.11497 | DR4 gives a *backward* check |
| Foustoul et al. 2025, CRTS+ZTF | 36 + 221 weaker | 2505.06656 | Itself the main retest |
| Huijse et al., Gaia DR3 | 181 → 13 (v2) | 2505.16884 | El-Badry, Hogg & Rix 2026, PASP (2509.10601) |

Methods warning: Lin, Charisi & Haiman (arXiv:2505.14778).

**Gap:** No one has tested pre-Gaia candidates against Gaia DR4 epoch photometry, because DR4 is not yet public. No pre-registered blind test of SMBHB ephemerides exists.

**Scoop risk:**
- The Gaia DPAC DR4 papers include a variable-AGN catalogue paper, but no listed SMBHB or periodicity paper.
- Huijse et al. plan a *blind DR4 periodicity search*, which is a different question.
- Differentiation: only published ephemerides are tested; the test is pre-registered; and Gaia is combined with CRTS and ZTF for a baseline of about 20 years. The RNAAS note should go out within about 2–3 weeks of release.

## 4. Hypotheses

- **H0:** For the population of candidates, the published sinusoid predicts DR4 photometry no better than a DRW model conditioned on the same pre-Gaia data. The median log Bayes factor is ≈ 0 or below.
- **H1:** A significant subset of candidates (more than expected under DRW simulations) are better predicted by their sinusoids.
- **Expectation:** a mostly negative result, consistent with El-Badry et al. 2026. Any survivors become high-priority SMBHB targets.

## 5. Data

| Source | Details |
|---|---|
| **Gaia DR4 `epoch_photometry`** | G/BP/RP per field-of-view transit; 2014-07-25 → 2020-01-20; ~40–100 epochs per quasar (estimate). Access: astroquery `Gaia.load_data(..., data_release='Gaia DR4', retrieval_type='EPOCH_PHOTOMETRY')` in batches of a few hundred; cross-match by TAP first. |
| Gaia DR4 `vari_agn`, `qso_candidates` | Context and cross-match |
| CRTS DR2 light curves | Pre-Gaia baseline (2005–2013) |
| ZTF DR24 (IRSA) | Post-Gaia baseline (2018 – Oct 2025) |
| Candidate tables | From the papers above (VizieR or paper tables): coordinates, period, phase/epoch, amplitude, band |

Gaia's limit is G ≈ 20.7, so faint PS1 and DES candidates drop out. The expected testable sample is about 150–300 candidates.

## 6. Method

**Before 2 Dec 2026 (pre-registration):**
1. Compile the candidate tables into `data/candidates.csv` with ephemeris parameters.
2. For each candidate, fit (a) the published sinusoid (+ DRW noise) and (b) a pure DRW to the pre-Gaia data, using celerite2 or EzTao.
3. Write the predictive distributions as functions of time for 2014–2020 under both models.
4. Commit `preregistration/` (model code, frozen parameters, scoring rule, sample cuts, and the conversion between photometric systems) and push to GitHub. **The commit timestamp precedes the data release.**

**After release:**
5. Cross-match to DR4 `source_id`s; download the epoch photometry.
6. Fit per-survey offsets and scale (Gaia G vs CRTS V vs ZTF g/r) using only pre-specified free parameters.
7. Score each candidate: log predictive likelihood of the DR4 points under each model, giving a Bayes factor.
8. Population test: compare the Bayes factor distribution with DRW simulations (simulated "fake candidates" selected the same way).
9. Secondary: a joint ~20-year CRTS + Gaia + ZTF fit for the survivors.

## 7. Timeline

| When | Tasks |
|---|---|
| **22 Nov – 1 Dec** | Compile candidates; download CRTS and ZTF data; build the models; **commit the pre-registration** |
| 2–5 Dec | DR4 release; cross-match; download epoch photometry (expect archive congestion) |
| 6–14 Dec | Score; population test; figures |
| 15–22 Dec | Write-up, mentor review, RNAAS submission, arXiv |

## 8. Venue

- **RNAAS:** fast, free, ideal for a timely result.
- A longer **MNRAS / ApJ** paper with the mentor, if many survivors or interesting methods.
- arXiv (mentor endorsement required for first-time submitters).

## 9. Where the mentor helps most

The SMBHB physics and the interpretation of survivors; Bayesian model comparison; arXiv endorsement.

## 10. Risks

| Risk | Mitigation |
|---|---|
| Scoop right after release | Pre-registration plus a narrow, specific question; fast RNAAS |
| Archive congestion | Prepare queries in advance; batch; retry |
| Too few Gaia epochs per source | Report which candidates are testable; use the Bayes factor, not a period search |
| Photometric-system mismatch | Pre-specified offset and scale parameters; sensitivity analysis |
| Epoch photometry not public for all sources | Check the DR4 documentation on release day; scale down if needed |

## 11. Key references

- Graham, M. J. et al. (2015). MNRAS 453, 1562. arXiv:1507.07603
- Charisi, M. et al. (2016). MNRAS 463, 2145. arXiv:1604.01020
- Liu, T. et al. (2019). ApJ 884, 36. arXiv:1906.08315
- Chen, Y.-C. et al. (2022/24). arXiv:2206.11497
- Foustoul et al. (2025). arXiv:2505.06656
- Huijse, P. et al. arXiv:2505.16884
- El-Badry, K., Hogg, D. W. & Rix, H.-W. (2026). PASP. arXiv:2509.10601
- Lin, Charisi & Haiman (2025). arXiv:2505.14778
- Vaughan, S. et al. (2016). MNRAS 461, 3145
- Kelly, B. C., Bechtold, J. & Siemiginowska, A. (2009). ApJ 698, 895 (DRW)
- Gaia Collaboration, DR4 documentation (2026)

> ⚠️ These references were compiled by an AI research assistant. Verify every citation before using it.
