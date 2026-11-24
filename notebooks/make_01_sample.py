"""Generate and execute notebooks/01_sample.ipynb (task T1.3).

Usage: .venv/bin/python notebooks/make_01_sample.py
"""
from pathlib import Path

import nbformat as nbf
from nbconvert.preprocessors import ExecutePreprocessor

HERE = Path(__file__).resolve().parent
cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s))  # noqa: E731
code = lambda s: cells.append(nbf.v4.new_code_cell(s))  # noqa: E731

md("""# 01 — Candidate sample, magnitude cuts and expected Gaia testability

This notebook summarises `data/candidates.csv` and the pre-registration inputs in
`results/frozen/`. **No Gaia epoch photometry is used**: the Gaia magnitude is the DR3 mean
`phot_g_mean_mag` (used only for the G < 20.7 cut), and Gaia cadences are *mock* cadences.
Power numbers are therefore planning estimates; the binding testability rule is evaluated on
the real DR4 timestamps after release (see `preregistration/prereg.md`).""")
code("""import sys, numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, '../src')
import config as C
pd.set_option('display.width', 160)
cand = pd.read_csv(C.CANDIDATES_CSV)
dr3 = pd.read_csv(C.RAW / 'dr3_mean_g.csv')
cand = cand.merge(dr3, on='cand_id', how='left')
fr = C.RESULTS / 'frozen'
summ = pd.read_csv(fr / 'discovery_summary.csv')
power = pd.read_csv(fr / 'power_mock.csv')
params = pd.read_csv(fr / 'frozen_params.csv')
len(cand)""")
md("## Candidates per sample, and how many pass G < 20.7 (DR3 mean G)")
code("""cand['G_lt_20p7'] = cand['dr3_g_mean'] < C.GAIA_G_LIMIT
cand['has_dr3'] = cand['dr3_status'] == 'ok'
t = cand.groupby('sample').agg(n=('cand_id','size'), final=('in_final_sample','sum'),
                               dr3_match=('has_dr3','sum'), G_lt_20p7=('G_lt_20p7','sum'))
t""")
code("""fin = cand[(cand['sample']!='Charisi2016') | cand['in_final_sample']]
fin.groupby('sample')['G_lt_20p7'].agg(['size','sum'])""")
md("## Why candidates drop out before scoring")
code("""summ.groupby(['sample','status']).size().unstack(fill_value=0).T""")
md("## Survey magnitude vs Gaia DR3 mean G")
code("""fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for s, g in cand.groupby('sample'):
    ax[0].scatter(g['mag'], g['dr3_g_mean'], s=10, label=s)
ax[0].plot([15, 22], [15, 22], 'k:', lw=1); ax[0].axhline(C.GAIA_G_LIMIT, color='r', lw=1)
ax[0].set_xlabel('published survey magnitude'); ax[0].set_ylabel('Gaia DR3 mean G'); ax[0].legend(fontsize=7)
ax[1].hist([g['dr3_g_mean'].dropna() for _, g in cand.groupby('sample')], bins=np.arange(14, 22, 0.25),
           stacked=True, label=list(cand.groupby('sample').groups))
ax[1].axvline(C.GAIA_G_LIMIT, color='r'); ax[1].set_xlabel('Gaia DR3 mean G'); ax[1].legend(fontsize=7)
plt.tight_layout()""")
md("## Expected power on the mock DR4 cadence")
code("""pw = power.merge(params[['cand_id','period_d','sin_semi_amp','D_sigma','test_win_lo','test_win_hi','dr3_g_mean']], on='cand_id')
pw['cycles_in_window'] = (pw.test_win_hi - pw.test_win_lo) / pw.period_d
pw.groupby('sample')[['power_p05','power_lnB3','null_frac_lnB3']].describe().T.round(3)""")
code("""rules = {'power(p<0.05) >= 0.5': pw.power_p05 >= 0.5,
         'power(p<0.05) >= 0.8': pw.power_p05 >= 0.8,
         'P(lnB>3 | S) >= 0.8': pw.power_lnB3 >= 0.8,
         'P(lnB>3 | S) >= 0.5': pw.power_lnB3 >= 0.5}
pd.DataFrame({k: pw[v].groupby('sample').size() for k, v in rules.items()}).fillna(0).astype(int)""")
code("""fig, ax = plt.subplots(1, 2, figsize=(10, 4))
sc = ax[0].scatter(pw.sin_semi_amp / pw.D_sigma, pw.power_p05, c=pw.cycles_in_window, s=14)
plt.colorbar(sc, ax=ax[0], label='cycles in test window')
ax[0].set_xlabel('sinusoid semi-amplitude / DRW sigma'); ax[0].set_ylabel('power (p<0.05)')
ax[1].scatter(pw.dr3_g_mean, pw.power_p05, s=14); ax[1].set_xlabel('DR3 mean G'); ax[1].set_ylabel('power (p<0.05)')
plt.tight_layout()""")
md("""## Notes
* Phases/T0 are not published by any source paper; Graham+2015 and Charisi+2016 publish no
  period errors. The frozen sinusoid phase and amplitude come from the discovery data only.
* The mock cadence (`src/gaia.py::mock_cadence`) is a planning assumption.""")

nb = nbf.v4.new_notebook(cells=cells, metadata={"kernelspec": {"name": "python3", "display_name": "Python 3"}})
ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb, {"metadata": {"path": str(HERE)}})
nbf.write(nb, HERE / "01_sample.ipynb")
print("wrote", HERE / "01_sample.ipynb")
