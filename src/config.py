"""Project-wide paths and constants.

Values that define the pre-registered analysis (priors, grids, thresholds) live
here so that the analysis code and preregistration/frozen_params.csv read the
same numbers. Do not change anything in the "FROZEN" block after the
pre-registration commit; add a dated amendment instead.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RAW = DATA / "raw"                # downloaded catalogue tables (git-ignored)
LC_DIR = DATA / "lightcurves"     # CRTS / PTF / ZTF light curves (git-ignored)
RESULTS = ROOT / "results"
PREREG = ROOT / "preregistration"
NOTES = ROOT / "notes"
CANDIDATES_CSV = DATA / "candidates.csv"   # tracked (exception in .gitignore)

for _d in (DATA, RAW, LC_DIR, RESULTS):
    _d.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# Gaia DR4 facts (ESA, checked 2026-11-24)
# --------------------------------------------------------------------------
GAIA_DR4_START_MJD = 56863.0   # 2014-07-25 (start of DR4 data, same as DR3)
GAIA_DR4_END_MJD = 58868.0     # 2020-01-20 (end of DR4 data)
GAIA_DR3_END_MJD = 57901.0     # 2017-05-28 (end of DR3 data; DR4-only interval starts here)
ZTF_START_MJD = 58194.0        # 2018-03-17, start of ZTF public survey
GAIA_G_LIMIT = 20.7            # practical faint limit used for testability

# --------------------------------------------------------------------------
# FROZEN analysis parameters (mirrored in preregistration/frozen_params.csv)
# --------------------------------------------------------------------------
# Hyperparameter fit to the discovery (pre-Gaia) light curve: MAP with priors
TAU_PRIOR_MEDIAN_D = 300.0     # DRW timescale prior median, observed frame [d]
TAU_PRIOR_SIGMA_LN = 1.0       # log-normal width
TAU_BOUNDS_D = (5.0, 1.0e4)
SIGMA_DRW_BOUNDS = (1e-3, 2.0) # DRW rms [mag]
JITTER_OLD_BOUNDS = (0.0, 0.5) # extra white noise added to discovery-survey errors [mag]

# Period marginalisation: Gauss-Hermite nodes over N(P_pub, sigma_P^2)
N_PERIOD_NODES = 7
SIGMA_P_FRAC_DEFAULT = 0.05    # fractional period error used when none is published

# Gaia nuisance parameters (identical in both models)
GAIA_OFFSET_PRIOR_SIGMA = 1.0  # [mag] Gaussian prior on Gaia-minus-survey zero point (mean 0)
SCALE_PRIOR_SIGMA_LN = 0.35    # log-normal prior on the G/survey amplitude ratio (median 1)
N_SCALE_NODES = 9
JITTER_GAIA_PRIOR_SIGMA = 0.05 # half-normal prior on extra Gaia white noise [mag]
N_JITTER_NODES = 6

# Gaia data handling
VISIT_GAP_D = 1.0              # transits closer than this are averaged into one visit
GAIA_OUTLIER_SIGMA = 5.0       # visits > 5 sigma (robust) from the window median are dropped

# Decision rules
LNB_THRESHOLD = 3.0            # individual "survivor" threshold on ln B
TESTABLE_POWER = 0.80          # fraction of injections with lnB > threshold for "testable"
P_SINGLE = 0.05                # per-candidate null p-value cut used in the H1 count
ALPHA_POP = 0.01               # population test significance
FDR_Q = 0.05                   # Benjamini-Hochberg level
N_NULL_SIMS = 1000             # accepted null simulations per candidate
N_INJ_SIMS = 300               # injection simulations per candidate
N_PHASE_SCRAMBLES = 200        # H4 control
RANDOM_SEED = 20261202
