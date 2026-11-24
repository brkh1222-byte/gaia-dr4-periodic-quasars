"""Unit tests for src/models.py: every fast formula is checked against brute-force
dense Gaussian algebra on small problems."""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import multivariate_normal

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import models as M  # noqa: E402


def _toy_lc(rng, n=60, two_surveys=False):
    t = np.sort(rng.uniform(53500, 56400, n))
    e = rng.uniform(0.03, 0.08, n)
    x = M.simulate_drw(t, 0.2, 250.0, rng)[0] + 18.0
    sv = (t > 55500).astype(int) if two_surveys else np.zeros(n, int)
    y = x + 0.3 * sv + e * rng.standard_normal(n)
    return M.LightCurve(t, y, e, sv)


def _dense_prior(lc, hyp, period, t_new=None):
    """Joint prior covariance of latent x at old (+new) times, incl. basis terms."""
    H = M.basis(lc.t, lc.survey, lc.n_surveys, period)
    lam = M.prior_cov_diag(lc.n_surveys, hyp.sigma_s if period else None)
    K = M.drw_cov(lc.t, lc.t, hyp.sigma, hyp.tau) + H @ np.diag(lam) @ H.T
    return K, H, lam


def test_factor_backends_agree():
    rng = np.random.default_rng(1)
    lc = _toy_lc(rng)
    f1 = M.Factor(lc.t, lc.e, 0.2, 300.0, 0.02, "celerite")
    f2 = M.Factor(lc.t, lc.e, 0.2, 300.0, 0.02, "dense")
    B = rng.standard_normal((len(lc.t), 3))
    assert np.allclose(f1.solve(B), f2.solve(B), rtol=1e-8, atol=1e-8)
    assert np.isclose(f1.logdet, f2.logdet, rtol=1e-10)


@pytest.mark.parametrize("two", [False, True])
@pytest.mark.parametrize("period", [None, 1200.0])
def test_log_evidence_matches_bruteforce(two, period):
    rng = np.random.default_rng(2)
    lc = _toy_lc(rng, two_surveys=two)
    hyp = M.Hyper(sigma=0.2, tau=300.0, jitter=0.02, sigma_s=0.15)
    le = M.log_evidence(lc, hyp.sigma, hyp.tau, hyp.jitter, period, hyp.sigma_s)
    K, _, _ = _dense_prior(lc, hyp, period)
    K += np.diag(lc.e**2 + hyp.jitter**2)
    ref = multivariate_normal(np.zeros(len(lc.y)), K).logpdf(lc.y)
    assert np.isclose(le, ref, rtol=1e-9, atol=1e-6)


@pytest.mark.parametrize("period", [None, 1100.0])
@pytest.mark.parametrize("two", [False, True])
def test_gaia_predictive_matches_bruteforce(period, two):
    """Fixed s, j and offset prior: compare with the conditional of a joint Gaussian."""
    rng = np.random.default_rng(3)
    lc = _toy_lc(rng, two_surveys=two)
    hyp = M.Hyper(sigma=0.2, tau=300.0, jitter=0.02, sigma_s=0.15)
    t_g = np.sort(rng.uniform(56900, 58800, 25))
    e_g = rng.uniform(0.01, 0.05, 25)
    s, j, sig_c = 1.3, 0.03, 0.7
    y_g = 18.0 + 0.1 * rng.standard_normal(25)

    # Brute force. Latent vector: beta (p), d_old (n), d_new (m), c (1).
    H = M.basis(lc.t, lc.survey, lc.n_surveys, period)
    lam = M.prior_cov_diag(lc.n_surveys, hyp.sigma_s if period else None)
    p, n, m = H.shape[1], len(lc.t), len(t_g)
    tt = np.concatenate([lc.t, t_g])
    Kd = M.drw_cov(tt, tt, hyp.sigma, hyp.tau)
    # y_old = H beta + d_old + eps ; y_g = beta_0 + c + s*(Hsin beta_sin + d_new) + eps_g
    Pg = np.zeros((m, p)); Pg[:, 0] = 1.0
    if period:
        w = 2 * np.pi / period
        Pg[:, lc.n_surveys:] = s * np.column_stack([np.cos(w * t_g), np.sin(w * t_g)])
    Bmat = np.vstack([H, Pg])                      # loading of beta
    Dmat = np.zeros((n + m, n + m)); Dmat[:n, :n] = np.eye(n); Dmat[n:, n:] = s * np.eye(m)
    cov = Bmat @ np.diag(lam) @ Bmat.T + Dmat @ Kd @ Dmat.T
    cov[n:, n:] += sig_c**2
    cov += np.diag(np.concatenate([lc.e**2 + hyp.jitter**2, e_g**2 + j**2]))
    full = multivariate_normal(np.zeros(n + m), cov).logpdf(np.concatenate([lc.y, y_g]))
    old = multivariate_normal(np.zeros(n), cov[:n, :n]).logpdf(lc.y)
    ref = full - old

    cd = M.condition(lc, hyp, t_g, period)
    got = M.gaia_logpdf_given(cd, y_g, e_g, lc.n_surveys, [s], [0.0], [j], [0.0], offset_sigma=sig_c)
    assert np.isclose(got, ref, rtol=1e-7, atol=1e-5)


def test_nuisance_weights_normalised():
    s, ls, j, lj = M.nuisance_nodes()
    assert np.isclose(np.exp(ls).sum(), 1) and np.isclose(np.exp(lj).sum(), 1)
    assert np.isclose(np.sum(np.exp(ls) * np.log(s)), 0, atol=1e-12)        # median 1
    assert j.min() == 0 and np.all(j >= 0)
    # half-normal second moment = sigma^2
    assert np.isclose(np.sum(np.exp(lj) * j**2), M.C.JITTER_GAIA_PRIOR_SIGMA**2, rtol=1e-6)


def test_simulate_drw_covariance():
    rng = np.random.default_rng(4)
    t = np.array([0.0, 50.0, 400.0, 30.0])
    x = M.simulate_drw(t, 0.3, 200.0, rng, n_draws=200000)
    emp = np.cov(x.T)
    assert np.allclose(emp, M.drw_cov(t, t, 0.3, 200.0), atol=3e-3)


def test_fit_recovers_drw_parameters():
    rng = np.random.default_rng(5)
    t = np.sort(rng.uniform(53500, 56400, 400))
    e = np.full(400, 0.03)
    y = 18 + M.simulate_drw(t, 0.25, 200.0, rng)[0] + e * rng.standard_normal(400)
    h = M.fit_hyper(M.LightCurve(t, y, e), "D")
    assert 0.15 < h.sigma < 0.4 and 80 < h.tau < 500


# ---------------------------------------------------------------------------
# score.py
# ---------------------------------------------------------------------------
import score as S  # noqa: E402


def _scorer_setup(seed=7, amp=0.25, n_g=30):
    rng = np.random.default_rng(seed)
    t = np.sort(rng.uniform(53500, 56400, 250))
    e = np.full(250, 0.04)
    P = 1300.0
    x = M.simulate_drw(t, 0.15, 200.0, rng)[0] + amp * np.sin(2 * np.pi * t / P + 0.7)
    lc = M.LightCurve(t, 18 + x + e * rng.standard_normal(250), e)
    hD = M.fit_hyper(lc, "D")
    hS = M.fit_hyper(lc, "S", period=P)
    t_g = np.sort(rng.uniform(57230, 58860, n_g))
    e_g = np.full(n_g, 0.03)
    return lc, hS, hD, t_g, e_g, rng


def test_scorer_matches_models_module():
    lc, hS, hD, t_g, e_g, rng = _scorer_setup()
    sc = S.CandidateScorer(lc, hS, hD, t_g, e_g)
    y = 18 + 0.1 * rng.standard_normal(len(t_g))
    ref = M.ln_bayes_factor(lc, hS, hD, t_g, y, e_g)
    assert np.isclose(sc.lnB(y[None])[0], ref, rtol=1e-8, atol=1e-6)


def test_null_and_injection_behave():
    """Under D, E[lnB] <= 0 (it is minus a KL divergence); under S, E[lnB] >= 0."""
    lc, hS, hD, t_g, e_g, rng = _scorer_setup()
    sc = S.CandidateScorer(lc, hS, hD, t_g, e_g)
    ln0 = sc.null_lnB(3000, rng)
    ln1 = sc.inj_lnB(3000, rng)
    assert ln0.mean() < 0 < ln1.mean()
    # null p-values of null draws are ~uniform
    ref = sc.null_lnB(3000, rng)
    p = np.array([S.null_pvalue(v, ref) for v in ln0[:500]])
    assert abs(np.mean(p < 0.05) - 0.05) < 0.03


def test_phase_rotation_zero_is_identity():
    lc, hS, hD, t_g, e_g, rng = _scorer_setup()
    a = S.CandidateScorer(lc, hS, hD, t_g, e_g)
    b = S.CandidateScorer(lc, hS, hD, t_g, e_g, phase_rotation=2 * np.pi)
    y = a.S.sample(5, rng)
    assert np.allclose(a.lnB(y), b.lnB(y), atol=1e-6)


def test_bh_and_binomial():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216])
    assert S.benjamini_hochberg(p, 0.05).sum() == 2
    assert np.isclose(S.binomial_excess_pvalue(0, 10, 0.05), 1.0)
