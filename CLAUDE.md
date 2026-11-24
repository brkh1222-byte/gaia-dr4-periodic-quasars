# CLAUDE.md

This is a research project by an independent astrophysics researcher: a pre-registered blind test of published periodic-quasar (supermassive black hole binary) candidates using Gaia DR4 epoch photometry (DR4 released 2 Dec 2026). Read PROPOSAL.md first. A mentor who works on black holes may co-author.

- **Deadline:** everything in `preregistration/` must be committed and pushed **before 2 Dec 2026**. After that, never edit those files; add dated amendments.
- **Blindness:** do not download or look at any Gaia DR4 data before the pre-registration commit exists.
- **Stack:** Python with astropy, astroquery (Gaia TAP and DataLink), celerite2 (DRW), numpy and pandas. Laptop-scale.
- **Provenance:** every candidate row must cite its source paper or table. Never invent ephemeris values.
- **Writing and integrity:** Claude drafts the paper; the researcher reviews and edits it; the paper carries an AI-use disclosure. Authorship is undecided: the draft uses the placeholders `[AUTHORS — TBD]` and `[AI-USE STATEMENT — TBD]`. Explain methods when implementing them. Verify every citation.
- **Environment:** `uv venv --python 3.11 .venv` then `uv pip install -r requirements.txt` (pinned). No conda.
