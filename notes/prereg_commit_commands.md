# Commands to run the moment the researcher approves the pre-registration

Deadline: before the Gaia DR4 release, 2 Dec 2026 ~12:00 CET (11:00 UTC). Do it at least a
day early. Nothing here has been run.

```sh
cd ~/Desktop/gaia-dr4-periodic-quasars

# 0. Sanity: tests pass, checksums match, nothing else is uncommitted
.venv/bin/python -m pytest -q tests
(cd preregistration && shasum -a 256 -c CHECKSUMS.txt)
git status --short          # only preregistration/* should be listed

# 1. (Choice 3, if confirmed) put the work on the default branch first
git checkout main && git merge --ff-only research-plan && git push origin main
#    if --ff-only refuses: git merge research-plan  (then push)

# 2. Commit and push the pre-registration
git add preregistration/
git commit -m "Pre-registration: frozen models, parameters, predictions and decision rules (before Gaia DR4)

Analysis code commit 9284e836a12399c934edb368783b747c8f73a9f9.

Co-Authored-By: Claude Fable 5.2 <noreply@anthropic.com>"
git push origin main            # or research-plan if not merged

# 3. Record the timestamped hash everywhere it is cited
git rev-parse HEAD              # -> PREREG_HASH; put it in paper/main.tex (\prereghash),
                                #    README.md and notes/STATUS.md, then commit and push
git tag -a prereg-v1 -m "Pre-registration, public before Gaia DR4" && git push origin prereg-v1

# 4. Optional but recommended: independent timestamps
#    - GitHub release from tag prereg-v1 (gh release create prereg-v1 --notes-file notes/prereg_summary.md)
#    - Zenodo DOI via the GitHub-Zenodo link (needs the researcher's Zenodo account)
#    - Software Heritage: https://archive.softwareheritage.org/save/ (paste the repo URL)
```

## If the researcher changes one of the three choices
- Testability 0.8 instead of 0.5: edit `TESTABLE_POWER_P05` in `src/analyze_dr4.py`, commit,
  then `.venv/bin/python src/make_prereg.py --code-hash $(git rev-parse HEAD) --power-rule 0.8`.
- Drop or pool the Graham PTF group: edit §2 and §6 of `preregistration/prereg_template.md`
  (pooling: run H1 with `--groups forecast,forecast_alt`), then rerun `make_prereg.py`.
- After any change: rerun step 0.

## Release day (after the pre-registration commit is public and DR4 is out)
```sh
.venv/bin/python src/gaia_xmatch.py --i-confirm-prereg-is-public          # table: gaiadr4.gaia_source
.venv/bin/python src/gaia_epochs.py --i-confirm-prereg-is-public          # 'Gaia DR4', EPOCH_PHOTOMETRY, DATAMODEL_STANDARD
.venv/bin/python src/analyze_dr4.py --groups forecast                     # primary
for s in dr4only unbinned pmarg widescale shuffle; do .venv/bin/python src/analyze_dr4.py --groups forecast --sensitivity $s; done
.venv/bin/python src/analyze_dr4.py --groups forecast_alt --tag dr4_alt
.venv/bin/python src/analyze_dr4.py --groups gap --tag dr4_gap
.venv/bin/python src/figures.py results/dr4
```
If DR4 table, release or column names differ, record the mapping in
`preregistration/amendment_2026-12-02.md` BEFORE scoring (the scripts take `--table`,
`--data-release`, `--retrieval-type`, `--data-structure`).
