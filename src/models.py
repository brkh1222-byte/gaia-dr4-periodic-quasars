"""Light-curve models for the pre-registered Gaia DR4 test (task T2.1).

Two models are compared. Both describe the latent magnitude x(t) of a quasar as

    x(t) = h(t)^T beta + d(t)

where d(t) is a damped random walk (DRW; Kelly et al. 2009) with covariance

    k(t, t') = sigma^2 exp(-|t - t'| / tau),

and beta are linear coefficients with Gaussian priors that are integrated out
analytically ("GP with explicit basis functions", Rasmussen & Williams 2006, §2.7):

* Model D (pure DRW):           h(t) = [1_survey(t)]              (one mean level per survey)
* Model S (sinusoid + DRW):     h(t) = [1_survey(t), cos(2 pi t/P), sin(2 pi t/P)]

In model S the sinusoid's amplitude and phase are the coefficients (a, b) with prior
N(0, sigma_s^2) each. Conditioning on the discovery light curve therefore gives the
posterior of the phase and amplitude at the published period P, which is exactly the
information a published period plus its discovery data carry. (No paper in our sample
publishes a phase or reference epoch; see notes/catalogue_provenance.md.) The period is
fixed at the published value in the primary analysis (researcher's decision, option A); a
sensitivity analysis integrates it over a grid of nodes with prior weights from the published
period and its error, reweighted by the evidence of the discovery data.

Hyperparameters (sigma, tau, extra white noise on the discovery data) are fitted to the
discovery data only, by maximum a posteriori with fixed priors (config.py), and frozen. The
sinusoid coefficients have a broad fixed prior (sigma_s = SIN_PRIOR_SIGMA = 1 mag), so their
posterior (phase and amplitude) is set by the discovery data.

Given frozen hyperparameters, `predict_gaia_logpdf` returns the log predictive density of
new (Gaia) magnitudes, y_G = mu0 + c + s * [x(t) - mu0] + noise, marginalised over the
nuisance parameters c (zero-point offset), s (amplitude ratio) and Gaia jitter j, which
are identical in both models.

Linear algebra: the discovery-data covariance is factorised once with celerite2 (O(N));
a dense backend is kept for testing.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import optimize
from scipy.linalg import cho_factor, cho_solve
from scipy.special import logsumexp

import celerite2
from celerite2 import terms

try:  # allow both `python -m src.models` and `sys.path.insert(0, "src")`
    from . import config as C
except ImportError:  # pragma: no cover
    import config as C

LN2PI = np.log(2.0 * np.pi)
MU_PRIOR_SIGMA = 10.0  # [mag] broad prior on per-survey mean levels (cancels in all comparisons)


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------
@dataclass
class LightCurve:
    """Discovery (pre-cutoff) light curve, possibly from several surveys."""

    t: np.ndarray            # MJD
    y: np.ndarray            # magnitude
    e: np.ndarray            # 1-sigma error [mag]
    survey: np.ndarray = None  # integer survey index per point (0 = reference survey)

    def __post_init__(self):
        self.t = np.asarray(self.t, float)
        self.y = np.asarray(self.y, float)
        self.e = np.asarray(self.e, float)
        if self.survey is None:
            self.survey = np.zeros(len(self.t), int)
        self.survey = np.asarray(self.survey, int)
        o = np.argsort(self.t, kind="stable")
        self.t, self.y, self.e, self.survey = self.t[o], self.y[o], self.e[o], self.survey[o]

    @property
    def n_surveys(self) -> int:
        return int(self.survey.max()) + 1

    @property
    def baseline(self) -> float:
        return float(self.t[-1] - self.t[0])


@dataclass
class Hyper:
    """Frozen hyperparameters of one model for one candidate."""

    sigma: float              # DRW rms [mag]
    tau: float                # DRW timescale [d], observed frame
    jitter: float             # extra white noise on discovery data [mag]
    sigma_s: float = 0.0      # sinusoid coefficient prior sd [mag] (model S only)
    period: float = np.nan    # published period [d] (model S only)
    sigma_p: float = np.nan   # period prior sd [d] (model S only)
    extra: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Basis functions and covariance pieces
# ---------------------------------------------------------------------------
def basis(t, survey, n_surveys, period=None, t_ref=0.0):
    """Design matrix H (N x p): survey indicator columns, then cos/sin if period given."""
    t = np.asarray(t, float)
    cols = [(np.asarray(survey) == k).astype(float) for k in range(n_surveys)]
    if period is not None:
        w = 2.0 * np.pi / period
        cols += [np.cos(w * (t - t_ref)), np.sin(w * (t - t_ref))]
    return np.column_stack(cols)


def prior_cov_diag(n_surveys, sigma_s=None):
    d = [MU_PRIOR_SIGMA**2] * n_surveys
    if sigma_s is not None:
        d += [sigma_s**2, sigma_s**2]
    return np.array(d)


def drw_cov(t1, t2, sigma, tau):
    return sigma**2 * np.exp(-np.abs(np.subtract.outer(t1, t2)) / tau)


class Factor:
    """Factorisation of K_y = K_drw(t,t) + diag(e^2 + jitter^2).

    backend='celerite' uses celerite2 (O(N)); backend='dense' uses a Cholesky.
    Exposes solve(B) = K_y^{-1} B and logdet.
    """

    def __init__(self, t, e, sigma, tau, jitter, backend="celerite"):
        self.t = np.asarray(t, float)
        diag = np.asarray(e, float) ** 2 + jitter**2
        self.backend = backend
        if backend == "celerite":
            gp = celerite2.GaussianProcess(terms.RealTerm(a=sigma**2, c=1.0 / tau), mean=0.0)
            gp.compute(self.t, diag=diag, quiet=False)
            self.gp = gp
            # log N(0 | 0, K) = -0.5 logdet - 0.5 N ln 2pi
            self.logdet = -2.0 * gp.log_likelihood(np.zeros_like(self.t)) - len(self.t) * LN2PI
        else:
            K = drw_cov(self.t, self.t, sigma, tau) + np.diag(diag)
            self.cf = cho_factor(K, lower=True)
            self.logdet = 2.0 * np.sum(np.log(np.diag(self.cf[0])))

    def solve(self, B):
        if self.backend == "celerite":
            return self.gp.apply_inverse(B)
        return cho_solve(self.cf, B)


# ---------------------------------------------------------------------------
# Evidence of the discovery data (used for hyperparameter fits and period weights)
# ---------------------------------------------------------------------------
def log_evidence(lc: LightCurve, sigma, tau, jitter, period=None, sigma_s=None,
                 backend="celerite", fac: Factor | None = None):
    """ln p(y | hyperparameters), with basis coefficients marginalised (Woodbury)."""
    if fac is None:
        fac = Factor(lc.t, lc.e, sigma, tau, jitter, backend)
    H = basis(lc.t, lc.survey, lc.n_surveys, period)
    lam = prior_cov_diag(lc.n_surveys, sigma_s if period is not None else None)
    KiH = fac.solve(H)
    Kiy = fac.solve(lc.y)
    A = np.diag(1.0 / lam) + H.T @ KiH          # p x p
    u = H.T @ Kiy
    cA = cho_factor(A, lower=True)
    quad = lc.y @ Kiy - u @ cho_solve(cA, u)
    logdet = fac.logdet + 2.0 * np.sum(np.log(np.diag(cA[0]))) + np.sum(np.log(lam))
    return -0.5 * (quad + logdet + len(lc.y) * LN2PI)


def log_prior_hyper(sigma, tau, jitter, sigma_s=None):
    """Weak MAP priors (config.py). Log-normal on tau; log-uniform on the rest."""
    lp = -0.5 * ((np.log(tau) - np.log(C.TAU_PRIOR_MEDIAN_D)) / C.TAU_PRIOR_SIGMA_LN) ** 2
    # half-normal on jitter (scale 0.1 mag) keeps it from absorbing the DRW
    lp += -0.5 * (jitter / 0.1) ** 2
    return lp


def fit_hyper(lc: LightCurve, model: str, period=None, sigma_p=None, n_starts=4,
              backend="celerite", rng=None, fit_sigma_s=False) -> Hyper:
    """MAP fit of hyperparameters to the discovery light curve.

    model 'D' fits (sigma, tau, jitter). Model 'S' fits (sigma, tau, jitter) with the period
    fixed at the published value and the sinusoid coefficients (a, b) given a broad, fixed
    prior N(0, SIN_PRIOR_SIGMA^2): the published sinusoid is taken as given and its phase and
    amplitude are estimated from the discovery data (researcher's option A). Fitting sigma_s
    by evidence maximisation (fit_sigma_s=True) is NOT used: it applies an Occam penalty that
    often removes the sinusoid altogether, i.e. it re-decides on the discovery data whether
    the sinusoid exists, which is the question the Gaia test is meant to answer.
    """
    fit_ss = (model == "S") and fit_sigma_s
    rng = np.random.default_rng(0) if rng is None else rng
    rms = float(np.std(lc.y))
    lo = [np.log(C.SIGMA_DRW_BOUNDS[0]), np.log(C.TAU_BOUNDS_D[0]), np.log(1e-4)]
    hi = [np.log(C.SIGMA_DRW_BOUNDS[1]), np.log(C.TAU_BOUNDS_D[1]), np.log(C.JITTER_OLD_BOUNDS[1])]
    if fit_ss:
        lo.append(np.log(1e-3)); hi.append(np.log(2.0))
    bounds = list(zip(lo, hi))

    def unpack(p):
        s, tau, j = np.exp(p[:3])
        ss = (np.exp(p[3]) if fit_ss else C.SIN_PRIOR_SIGMA) if model == "S" else None
        return s, tau, j, ss

    def nlp(p):
        s, tau, j, ss = unpack(p)
        try:
            le = log_evidence(lc, s, tau, j, period if model == "S" else None, ss, backend)
        except Exception:
            return 1e25
        val = -(le + log_prior_hyper(s, tau, j, ss))
        return val if np.isfinite(val) else 1e25

    best = None
    for k in range(n_starts):
        p0 = [np.log(max(rms, 0.02)), np.log(C.TAU_PRIOR_MEDIAN_D), np.log(0.02)]
        if fit_ss:
            p0.append(np.log(max(rms, 0.02)))
        p0 = np.array(p0) + (0 if k == 0 else rng.normal(0, 0.7, len(p0)))
        p0 = np.clip(p0, np.array(lo) + 1e-6, np.array(hi) - 1e-6)
        r = optimize.minimize(nlp, p0, method="L-BFGS-B", bounds=bounds)
        if best is None or r.fun < best.fun:
            best = r
    s, tau, j, ss = unpack(best.x)
    h = Hyper(sigma=s, tau=tau, jitter=j)
    if model == "S":
        h.sigma_s, h.period = ss, float(period)
        h.sigma_p = float(sigma_p) if sigma_p is not None else C.SIGMA_P_FRAC_DEFAULT * period
    h.extra = {"neg_log_post": float(best.fun), "converged": bool(best.success)}
    return h


# ---------------------------------------------------------------------------
# Period nodes for model S
# ---------------------------------------------------------------------------
def period_nodes(period, sigma_p, n=None):
    """Uniform grid over +-3 sigma_p with Gaussian prior weights (normalised)."""
    n = C.N_PERIOD_NODES if n is None else n
    x = np.linspace(-3.0, 3.0, n)
    P = period + sigma_p * x
    P = P[P > 0]
    x = (P - period) / sigma_p
    logw = -0.5 * x**2
    return P, logw - logsumexp(logw)


# ---------------------------------------------------------------------------
# Gaia nuisance quadrature
# ---------------------------------------------------------------------------
def nuisance_nodes():
    """Quadrature nodes and log-weights for the amplitude ratio s and the Gaia jitter j."""
    xs, ws = np.polynomial.hermite_e.hermegauss(C.N_SCALE_NODES)
    s_nodes = np.exp(C.SCALE_PRIOR_SIGMA_LN * xs)
    s_logw = np.log(ws / ws.sum())
    # half-normal via the non-negative Gauss-Hermite nodes (odd order includes 0)
    m = 2 * C.N_JITTER_NODES - 1
    xj, wj = np.polynomial.hermite_e.hermegauss(m)
    keep = xj >= -1e-12
    wj = np.where(np.abs(xj[keep]) < 1e-12, wj[keep], 2 * wj[keep])
    j_nodes = C.JITTER_GAIA_PRIOR_SIGMA * xj[keep]
    j_logw = np.log(wj / wj.sum())
    return s_nodes, s_logw, j_nodes, j_logw


# ---------------------------------------------------------------------------
# Predictive distribution of the latent light curve at new times
# ---------------------------------------------------------------------------
@dataclass
class Conditioned:
    """Pieces of the posterior predictive needed for Gaia scoring (one period node)."""

    A_y: np.ndarray     # K_no K_y^{-1} y_o                (n_new,)
    AH: np.ndarray      # K_no K_y^{-1} H_o                (n_new, p)
    C0: np.ndarray      # K_nn - K_no K_y^{-1} K_on        (n_new, n_new)
    m: np.ndarray       # posterior mean of beta           (p,)
    S: np.ndarray       # posterior covariance of beta     (p, p)
    Hn_sin: np.ndarray  # sinusoid basis at new times      (n_new, 2) or (n_new, 0)
    ref_cols: np.ndarray  # indices of the reference-survey mean column
    log_ev: float       # ln p(y_o | hyper, node)


def condition(lc: LightCurve, hyp: Hyper, t_new, period=None, fac: Factor | None = None,
              backend="celerite", KinvKon=None) -> Conditioned:
    """Condition the model on the discovery data and return predictive pieces at t_new.

    Gaia magnitudes are linked to the reference survey (index 0).
    """
    t_new = np.asarray(t_new, float)
    if fac is None:
        fac = Factor(lc.t, lc.e, hyp.sigma, hyp.tau, hyp.jitter, backend)
    use_sin = period is not None
    H = basis(lc.t, lc.survey, lc.n_surveys, period)
    lam = prior_cov_diag(lc.n_surveys, hyp.sigma_s if use_sin else None)
    KiH = fac.solve(H)
    Kiy = fac.solve(lc.y)
    Sinv = np.diag(1.0 / lam) + H.T @ KiH
    cS = cho_factor(Sinv, lower=True)
    S = cho_solve(cS, np.eye(len(lam)))
    m = S @ (H.T @ Kiy)
    # evidence (same algebra as log_evidence)
    u = H.T @ Kiy
    quad = lc.y @ Kiy - u @ cho_solve(cS, u)
    logdet = fac.logdet + 2.0 * np.sum(np.log(np.diag(cS[0]))) + np.sum(np.log(lam))
    log_ev = -0.5 * (quad + logdet + len(lc.y) * LN2PI)

    Kno = drw_cov(t_new, lc.t, hyp.sigma, hyp.tau)
    if KinvKon is None:
        KinvKon = fac.solve(Kno.T)            # (n_old, n_new)
    A_y = KinvKon.T @ lc.y
    AH = KinvKon.T @ H
    C0 = drw_cov(t_new, t_new, hyp.sigma, hyp.tau) - Kno @ KinvKon
    C0 = 0.5 * (C0 + C0.T)
    if use_sin:
        w = 2.0 * np.pi / period
        Hn_sin = np.column_stack([np.cos(w * t_new), np.sin(w * t_new)])
    else:
        Hn_sin = np.zeros((len(t_new), 0))
    return Conditioned(A_y, AH, C0, m, S, Hn_sin, np.array([0]), log_ev)


def gaia_moments(cd: Conditioned, s: float, n_surveys: int):
    """Mean and covariance of y_G (before adding the offset and noise terms).

    y_G = mu_0 + c + s * [sin_part + d_new] + eps  with d_new = A y_o - A H beta + w.
    """
    n_new = cd.C0.shape[0]
    p = len(cd.m)
    P = np.zeros((n_new, p))
    P[:, 0] = 1.0                                  # reference-survey mean level, unscaled
    if cd.Hn_sin.shape[1]:
        P[:, n_surveys:n_surveys + 2] = s * cd.Hn_sin
    Q = P - s * cd.AH
    mean = s * cd.A_y + Q @ cd.m
    cov = Q @ cd.S @ Q.T + s**2 * cd.C0
    return mean, cov


def gaia_logpdf_given(cd: Conditioned, y_g, e_g, n_surveys, s_nodes, s_logw, j_nodes, j_logw,
                      offset_sigma=None):
    """ln p(y_G | discovery data) marginalised over offset c (analytic), s and j (quadrature)."""
    offset_sigma = C.GAIA_OFFSET_PRIOR_SIGMA if offset_sigma is None else offset_sigma
    y_g = np.asarray(y_g, float)
    e2 = np.asarray(e_g, float) ** 2
    n = len(y_g)
    one = np.ones(n)
    out = []
    for s, ls in zip(s_nodes, s_logw):
        mean, cov = gaia_moments(cd, s, n_surveys)
        cov = cov + offset_sigma**2 * np.outer(one, one)
        r = y_g - mean
        for j, lj in zip(j_nodes, j_logw):
            K = cov + np.diag(e2 + j**2)
            cf = cho_factor(K, lower=True)
            ll = -0.5 * (r @ cho_solve(cf, r) + 2 * np.sum(np.log(np.diag(cf[0]))) + n * LN2PI)
            out.append(ls + lj + ll)
    return logsumexp(out)


# ---------------------------------------------------------------------------
# Top-level scoring helpers
# ---------------------------------------------------------------------------
def predictive_logpdf_D(lc, hypD, t_g, y_g, e_g, backend="celerite"):
    nodes = nuisance_nodes()
    cd = condition(lc, hypD, t_g, None, backend=backend)
    return gaia_logpdf_given(cd, y_g, e_g, lc.n_surveys, *nodes)


def predictive_logpdf_S(lc, hypS, t_g, y_g, e_g, backend="celerite", return_weights=False,
                        marginalise_period=False):
    """Model S predictive. Primary analysis: period fixed at the published value.

    With marginalise_period=True (sensitivity analysis): mixture over period nodes with
    weights proportional to prior x evidence of the discovery data.
    """
    nodes = nuisance_nodes()
    fac = Factor(lc.t, lc.e, hypS.sigma, hypS.tau, hypS.jitter, backend)
    Kno = drw_cov(np.asarray(t_g, float), lc.t, hypS.sigma, hypS.tau)
    KinvKon = fac.solve(Kno.T)
    if marginalise_period:
        Ps, logw = period_nodes(hypS.period, hypS.sigma_p)
    else:
        Ps, logw = np.array([hypS.period]), np.array([0.0])
    lls, lev = [], []
    for P in Ps:
        cd = condition(lc, hypS, t_g, P, fac=fac, KinvKon=KinvKon)
        lev.append(cd.log_ev)
        lls.append(gaia_logpdf_given(cd, y_g, e_g, lc.n_surveys, *nodes))
    post = logw + np.array(lev)
    post -= logsumexp(post)
    val = logsumexp(post + np.array(lls))
    if return_weights:
        return val, Ps, np.exp(post)
    return val


def ln_bayes_factor(lc, hypS, hypD, t_g, y_g, e_g, backend="celerite"):
    """ln B = ln p(D_new | D_old, S) - ln p(D_new | D_old, D)."""
    return (predictive_logpdf_S(lc, hypS, t_g, y_g, e_g, backend)
            - predictive_logpdf_D(lc, hypD, t_g, y_g, e_g, backend))


# ---------------------------------------------------------------------------
# Simulation helpers (used by tests and the null simulation)
# ---------------------------------------------------------------------------
def simulate_drw(t, sigma, tau, rng, n_draws=1):
    """Exact DRW (Ornstein-Uhlenbeck) draws at sorted or unsorted times; shape (n_draws, N)."""
    t = np.asarray(t, float)
    o = np.argsort(t, kind="stable")
    ts = t[o]
    dt = np.diff(ts)
    a = np.exp(-dt / tau)
    x = np.empty((n_draws, len(ts)))
    x[:, 0] = sigma * rng.standard_normal(n_draws)
    z = rng.standard_normal((n_draws, len(ts) - 1))
    sd = sigma * np.sqrt(1.0 - a**2)
    for i in range(1, len(ts)):
        x[:, i] = a[i - 1] * x[:, i - 1] + sd[i - 1] * z[:, i - 1]
    out = np.empty_like(x)
    out[:, o] = x
    return out


def predictive_moments(lc, hyp, t_new, period=None, backend="celerite"):
    """Mean and sd of the latent light curve x(t_new) in the reference survey (for plots)."""
    cd = condition(lc, hyp, t_new, period, backend=backend)
    mean, cov = gaia_moments(cd, 1.0, lc.n_surveys)
    return mean, np.sqrt(np.clip(np.diag(cov), 0, None))
