# CLAUDE.md

This is a research project by a 17-year-old astrophysics student: a pre-registered blind test of published periodic-quasar (supermassive black hole binary) candidates using Gaia DR4 epoch photometry (DR4 released 2 Dec 2026). Read PROPOSAL.md first. The student's mentor works on black holes and may co-author.

- **Deadline:** everything in `preregistration/` must be committed and pushed **before 2 Dec 2026**. After that, never edit those files; add dated amendments.
- **Blindness:** do not download or look at any Gaia DR4 data before the pre-registration commit exists.
- **Stack:** Python with astropy, astroquery (Gaia TAP and DataLink), celerite2 (DRW), numpy and pandas. Laptop-scale.
- **Provenance:** every candidate row must cite its source paper or table. Never invent ephemeris values.
- **Integrity:** the student must understand and write the paper themselves. Explain methods when implementing them. Verify every citation.
