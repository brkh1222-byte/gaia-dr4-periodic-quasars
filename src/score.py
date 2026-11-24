"""Scoring of Gaia epochs against the two frozen models (task T2.2).

For one candidate, `CandidateScorer` holds both models conditioned on the discovery light
curve (frozen hyperparameters) and the Gaia test epochs (times and errors). It returns

  ln B = ln p(y_G | D_old, S) - ln p(y_G | D_old, D)

for the observed Gaia magnitudes, and can also *simulate* Gaia magnitudes from either
model's predictive distribution. Simulating from model D gives the conditional null
distribution of ln B (the p-value of the candidate); simulating from model S gives the
injection distribution (the power, which defines testability).

The free nuisance parameters are the same in both models and fixed in advance:
  c  Gaia-minus-survey zero point   ~ N(0, 1.0^2) mag   (analytic marginalisation)
  s  G / survey amplitude ratio     ~ LogNormal(0, 0.35) (Gauss-Hermite quadrature)
  j  extra Gaia white noise         ~ HalfNormal(0.05) mag (Gauss-Hermite quadrature)

Primary analysis (researcher's decision, option A): the period is fixed at the published
value; phase and amplitude are the posterior of the sinusoid coefficients given the
discovery data only. A period-marginalised variant exists for the sensitivity analysis.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.special import logsumexp

try:
    from . import config as C
    from . import models as M
except ImportError:  # pragma: no cover
    import config as C
    import models as M

LN2PI = np.log(2 * np.pi)


class _Component:
    """Gaussian predictive y_G ~ N(mean_k, cov_k) for one (period node, s, j) combination."""

    __slots__ = ("mean", "L", "logdet", "logw", "cov")

    def __init__(self, mean, cov, logw):
        self.mean, self.cov, self.logw = mean, cov, logw
        self.L = np.linalg.cholesky(cov)
        self.logdet = 2.0 * np.sum(np.log(np.diag(self.L)))

    def logpdf(self, Y):
        """Y: (n_draws, n) -> (n_draws,)"""
        R = (Y - self.mean).T                       # (n, n_draws)
        Z = solve_triangular(self.L, R, lower=True)
        return -0.5 * (np.sum(Z * Z, axis=0) + self.logdet + len(self.mean) * LN2PI)


class ModelPredictive:
    """Mixture-of-Gaussians predictive of one model at the Gaia epochs."""

    def __init__(self, lc, hyp, t_g, e_g, periods=None, period_logw=None,
                 phase_rotation=0.0, offset_sigma=None, backend="celerite", scale_sigma_ln=None):
        offset_sigma = C.GAIA_OFFSET_PRIOR_SIGMA if offset_sigma is None else offset_sigma
        t_g = np.asarray(t_g, float)
        e2 = np.asarray(e_g, float) ** 2
        n = len(t_g)
        s_nodes, s_logw, j_nodes, j_logw = M.nuisance_nodes(scale_sigma_ln)
        fac = M.Factor(lc.t, lc.e, hyp.sigma, hyp.tau, hyp.jitter, backend)
        Kno = M.drw_cov(t_g, lc.t, hyp.sigma, hyp.tau)
        KinvKon = fac.solve(Kno.T)
        if periods is None:
            nodes = [(None, 0.0)]
        else:
            nodes = list(zip(periods, period_logw))
        conds, lev = [], []
        for P, lw in nodes:
            cd = M.condition(lc, hyp, t_g, P, fac=fac, KinvKon=KinvKon)
            if P is not None and phase_rotation:
                cd = _rotate_phase(cd, lc.n_surveys, phase_rotation)
            conds.append((cd, lw))
            lev.append(cd.log_ev)
        # posterior weights of period nodes: prior x evidence of discovery data
        lw = np.array([w for _, w in conds]) + np.array(lev)
        lw -= logsumexp(lw)
        self.period_weights = np.exp(lw)
        self.periods = [P for P, _ in nodes]
        self.comps = []
        one = np.ones((n, n))
        for (cd, _), lwp in zip(conds, lw):
            for s, ls in zip(s_nodes, s_logw):
                mean, cov = M.gaia_moments(cd, s, lc.n_surveys)
                cov = cov + offset_sigma**2 * one
                for j, lj in zip(j_nodes, j_logw):
                    K = cov + np.diag(e2 + j**2)
                    self.comps.append(_Component(mean, K, lwp + ls + lj))
        self.logw = np.array([c.logw for c in self.comps])
        self.n = n

    def logpdf(self, Y):
        Y = np.atleast_2d(Y)
        L = np.stack([c.logpdf(Y) for c in self.comps], axis=0)   # (n_comp, n_draws)
        return logsumexp(L + self.logw[:, None], axis=0)

    def sample(self, n_draws, rng):
        """Draw Gaia light curves from the predictive (nuisances drawn from their priors)."""
        w = np.exp(self.logw - logsumexp(self.logw))
        k = rng.choice(len(self.comps), size=n_draws, p=w)
        out = np.empty((n_draws, self.n))
        for ki in np.unique(k):
            idx = np.where(k == ki)[0]
            c = self.comps[ki]
            z = rng.standard_normal((self.n, len(idx)))
            out[idx] = (c.mean[:, None] + c.L @ z).T
        return out


def _rotate_phase(cd, n_surveys, phi):
    """Shift the predicted sinusoid by phase phi at the Gaia epochs (H4 control).

    cos(wt+phi) = cos wt cos phi - sin wt sin phi ; sin(wt+phi) = sin wt cos phi + cos wt sin phi.
    Only the forecast is shifted; the fit to the discovery data is unchanged.
    """
    R = np.array([[np.cos(phi), np.sin(phi)], [-np.sin(phi), np.cos(phi)]])
    return M.Conditioned(cd.A_y, cd.AH, cd.C0, cd.m, cd.S, cd.Hn_sin @ R,
                         cd.ref_cols, cd.log_ev)


class CandidateScorer:
    """Both models for one candidate at its Gaia test epochs."""

    def __init__(self, lc, hypS, hypD, t_g, e_g, period=None, sigma_p=None,
                 marginalise_period=False, phase_rotation=0.0, backend="celerite",
                 scale_sigma_ln=None):
        period = hypS.period if period is None else period
        if marginalise_period:
            sp = hypS.sigma_p if sigma_p is None else sigma_p
            Ps, lw = M.period_nodes(period, sp)
        else:
            Ps, lw = np.array([period]), np.array([0.0])
        self.S = ModelPredictive(lc, hypS, t_g, e_g, Ps, lw, phase_rotation, backend=backend,
                                 scale_sigma_ln=scale_sigma_ln)
        self.D = ModelPredictive(lc, hypD, t_g, e_g, None, None, backend=backend,
                                 scale_sigma_ln=scale_sigma_ln)

    def lnB(self, Y):
        return self.S.logpdf(Y) - self.D.logpdf(Y)

    def null_lnB(self, n_draws, rng):
        """ln B for Gaia data simulated from model D (conditional null)."""
        return self.lnB(self.D.sample(n_draws, rng))

    def inj_lnB(self, n_draws, rng):
        """ln B for Gaia data simulated from model S (power / injection)."""
        return self.lnB(self.S.sample(n_draws, rng))


def null_pvalue(lnb_obs, lnb_null):
    """One-sided p = (1 + #{null >= obs}) / (1 + N)."""
    lnb_null = np.asarray(lnb_null)
    return (1 + np.sum(lnb_null >= lnb_obs)) / (1 + len(lnb_null))


def benjamini_hochberg(p, q=None):
    """Boolean mask of discoveries at FDR level q."""
    q = C.FDR_Q if q is None else q
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    thresh = q * (np.arange(1, n + 1) / n)
    below = p[o] <= thresh
    k = np.max(np.where(below)[0]) + 1 if below.any() else 0
    mask = np.zeros(n, bool)
    mask[o[:k]] = True
    return mask


def binomial_excess_pvalue(k, n, rate):
    """One-sided P(X >= k) for X ~ Binomial(n, rate)."""
    from scipy.stats import binom
    return float(binom.sf(k - 1, n, rate))


def fraction_posterior(lnb_obs, null_draws, inj_draws, f_grid=None, bw=None):
    """Posterior on the fraction f of true periodic sources (flat prior).

    L(f) = prod_i [ f * p_inj,i(lnB_i) + (1-f) * p_null,i(lnB_i) ], with p_* estimated by a
    Gaussian KDE of each candidate's simulated ln B (Silverman bandwidth unless given).
    """
    from scipy.stats import gaussian_kde
    f_grid = np.linspace(0, 1, 1001) if f_grid is None else f_grid
    logL = np.zeros_like(f_grid)
    for x, dn, di in zip(lnb_obs, null_draws, inj_draws):
        pn = gaussian_kde(dn, bw_method=bw)(x)[0] + 1e-300
        pi = gaussian_kde(di, bw_method=bw)(x)[0] + 1e-300
        logL += np.log(f_grid * pi + (1 - f_grid) * pn)
    post = np.exp(logL - logL.max())
    post /= np.trapezoid(post, f_grid)
    cdf = np.cumsum(post) * (f_grid[1] - f_grid[0])
    return dict(f=f_grid, post=post, median=float(np.interp(0.5, cdf, f_grid)),
                upper95=float(np.interp(0.95, cdf, f_grid)))
