# Novelty and facts check — 2026-11-24

Web-based check (literature only; no Gaia epoch photometry was opened). arXiv listing
pages rate-limited us (HTTP 429), and ADS/Google Scholar were not reached, so the
sweep is **not exhaustive**. Repeat on 2026-11-30 before the pre-registration commit.

## Gaia release facts (ESA pages, checked 2026-11-24)
- DR4 release: **2 December 2026, ~12:00 CET** (banner on
  https://www.cosmos.esa.int/web/gaia/data-release-4; https://www.cosmos.esa.int/web/gaia/release).
- DR4 data span: 25 Jul 2014 10:30 UTC – 20 Jan 2020 22:00 UTC (66 months).
  `epoch_photometry` ~1.22e11 rows (7.2 TB); `vari_agn` 1,966,371 rows.
  https://www.cosmos.esa.int/web/gaia/dr4
- DR3 data span: 25 Jul 2014 10:30 UTC – 28 May 2017 08:44 UTC (34 months).
  DR3 `vari_agn`: **872,228** AGN with public G/BP/RP epoch photometry (mean 39 transits
  over 925 d). https://www.cosmos.esa.int/web/gaia/dr3 ; arXiv:2206.05681.
  -> Some candidates may already have public DR3 light curves for 2014.6–2017.4.
  The pre-registration forbids inspecting them and adds a DR4-only (after 2017-05-28)
  robustness test.

## Is anyone already doing this?
- **No paper, preprint or announcement found** that tests published CRTS/PTF/PS1/ZTF/
  Foustoul ephemerides against Gaia DR3 or DR4 photometry, and none that pre-registers
  predictions for a periodic-quasar sample.
- ESA DR4 papers list (95 titles, https://www.cosmos.esa.int/web/gaia/dr4-papers): no
  periodic-AGN / SMBHB paper. Closest: "Gaia DR4: The second Gaia catalogue of variable
  active galactic nuclei" (content not visible; may discuss periodicity); "Time delays in
  lensed quasars"; arXiv:2605.23388 (extragalactic classification performance).
- Huijse et al. 2026, A&A 715, A92 (arXiv:2505.16884; retitled "A search for periodic AGN
  variability in Gaia DR3"): 13/377,128 pass, consistent with false positives; says the
  pipeline is "readily applicable" to DR4 -> most likely fast DR4 periodicity search. It
  searches for new periods; it does not test published ephemerides.
- El-Badry, Hogg & Rix 2026 (arXiv:2509.10601): Gaia DR3 sinusoids -> ZTF; opposite
  direction, Gaia-selected sample.
- **Unverified lead:** a search snippet suggests a DR3 periodicity preprint cites a
  "Gaia Collaboration (2026)" paper describing candidates found without red-noise
  modelling as "overwhelmingly false positives" (possibly the DR4 variable-AGN paper).
  Not confirmed (aanda.org returned 403). Follow up.
- Not overlapping, noted for context: Agarwal et al. (NANOGrav) arXiv:2508.16534 (targeted
  PTA search, 114 candidates); GRAVITY Collab. 2026 arXiv:2610.14174; arXiv:2606.05343
  (PKS 2131-021 eccentricity); arXiv:2504.04278 (PKS 2131-021 optical red-noise null);
  arXiv:2603.09275 (4C 50.43, red noise).

## Conclusion
The gap stands as described in docs/research-plan.md §A.4: the test of fixed published
ephemerides against Gaia, pre-registered, appears unclaimed. The main scoop risk remains a
fast DR4 periodicity search (Huijse group or El-Badry/Hogg/Rix), which would answer a
different question.

## Other facts gathered
- RNAAS: 1,500 words including title, headers, captions and references (abstract ≤150
  words), one figure *or* one table; LaTeX preferably AASTeX v7+.
  https://journals.aas.org/research-note-preparation-guidelines/
- CRTS: old DR2 cone search (nunuku.caltech.edu) returns 404. Catalina Surveys portal
  http://nesssi.cacr.caltech.edu/DataRelease/ responds (single-position forms only; no
  bulk interface).
