"""Two-sample MR methods re-implemented in Python to mirror TwoSampleMR / MR-PRESSO (R).
Validated against the original R logs (see validate_against_R.py)."""
import numpy as np
from scipy import stats

def _wls(x, y, w):
    X = np.column_stack([np.ones_like(x), x]); W = np.diag(w)
    XtWX = X.T @ W @ X; inv = np.linalg.inv(XtWX)
    beta = inv @ X.T @ W @ y
    res = y - X @ beta; n = len(y)
    sigma2 = np.sum(w * res**2) / (n - 2)
    se = np.sqrt(np.diag(inv) * sigma2)
    return beta, se, np.sqrt(sigma2), res

def ivw(bx, by, sey):
    n = len(bx); w = 1 / sey**2
    S = np.sum(w * bx**2); b = np.sum(w * bx * by) / S
    Q = np.sum(w * (by - b * bx)**2); sigma = np.sqrt(Q / (n - 1))
    se = np.sqrt(1 / S) * max(1.0, sigma)             # TwoSampleMR: se_lm / min(1, sigma)
    p = 2 * stats.norm.sf(abs(b / se))              # matches the p-value reported by the R log (normal)
    return dict(b=b, se=se, pval=p, nsnp=n, Q=Q, Q_df=n - 1, Q_pval=stats.chi2.sf(Q, n - 1))

def ivw_fixed(bx, by, sey):
    w = 1 / sey**2; S = np.sum(w * bx**2); b = np.sum(w * bx * by) / S
    return b, np.sqrt(1 / S)

def egger(bx, by, sey):
    n = len(bx); s = np.sign(bx); s[s == 0] = 1
    x = np.abs(bx); y = by * s; w = 1 / sey**2
    beta, se, sigma, res = _wls(x, y, w)
    slope_se = se[1] / min(1.0, sigma); int_se_trunc = se[0] / min(1.0, sigma)
    Q = sigma**2 * (n - 2)
    out = dict(b=beta[1], se=slope_se, pval=2 * stats.t.sf(abs(beta[1] / slope_se), n - 2), nsnp=n,
               Q=Q, Q_df=n - 2, Q_pval=stats.chi2.sf(Q, n - 2))
    # TwoSampleMR::mr_pleiotropy_test reports the raw lm standard error for the intercept
    out.update(intercept=beta[0], intercept_se=se[0], intercept_pval=2 * stats.t.sf(abs(beta[0] / se[0]), n - 2),
               intercept_se_trunc=int_se_trunc)
    return out

def _wmedian(b_iv, w):
    o = np.argsort(b_iv); b = b_iv[o]; w = w[o]
    ws = (np.cumsum(w) - 0.5 * w) / np.sum(w)
    below = np.max(np.where(ws < 0.5)[0])
    return b[below] + (b[below + 1] - b[below]) * (0.5 - ws[below]) / (ws[below + 1] - ws[below])

def weighted_median(bx, by, sex, sey, nboot=1000, rng=None):
    rng = rng or np.random.default_rng(1)
    b_iv = by / bx
    VBj = sey**2 / bx**2 + (by**2 * sex**2) / bx**4
    w = 1 / VBj
    b = _wmedian(b_iv, w)
    meds = np.empty(nboot)
    for i in range(nboot):
        bxb = rng.normal(bx, sex); byb = rng.normal(by, sey)
        meds[i] = _wmedian(byb / bxb, w)
    se = np.std(meds, ddof=1)
    return dict(b=b, se=se, pval=2 * stats.norm.sf(abs(b / se)), nsnp=len(bx))

def _mad(x):
    return 1.4826 * np.median(np.abs(x - np.median(x)))

def _mode_est(b_iv, se_iv, weighted, phi=1.0):
    n = len(b_iv)
    s = 0.9 * min(np.std(b_iv, ddof=1), _mad(b_iv)) / n**0.2
    w = (se_iv**-2 / np.sum(se_iv**-2)) if weighted else np.full(n, 1.0 / n)
    h = max(1e-8, s * phi)
    grid = np.linspace(b_iv.min() - 3 * h, b_iv.max() + 3 * h, 512)
    dens = (w[None, :] * np.exp(-0.5 * ((grid[:, None] - b_iv[None, :]) / h)**2)).sum(axis=1)
    return grid[np.argmax(dens)]

def mode(bx, by, sex, sey, weighted, nboot=1000, phi=1.0, rng=None):
    rng = rng or np.random.default_rng(2)
    b_iv = by / bx; se_iv = sey / np.abs(bx)
    b = _mode_est(b_iv, se_iv, weighted, phi)
    boots = np.empty(nboot)
    for i in range(nboot):
        bxb = rng.normal(bx, sex); byb = rng.normal(by, sey)
        ivb = byb / bxb; seb = np.sqrt(sey**2 / bxb**2 + (byb**2 / bxb**4) * sex**2)
        boots[i] = _mode_est(ivb, seb, weighted, phi)
    se = _mad(boots)
    return dict(b=b, se=se, pval=2 * stats.t.sf(abs(b / se), len(bx) - 1), nsnp=len(bx))

def run_all(bx, by, sex, sey, nboot=1000, seed=1):
    rng = np.random.default_rng(seed)
    r = {}
    r['MR Egger'] = egger(bx, by, sey)
    r['Weighted median'] = weighted_median(bx, by, sex, sey, nboot, rng)
    r['Inverse variance weighted'] = ivw(bx, by, sey)
    r['Simple mode'] = mode(bx, by, sex, sey, False, nboot, rng=rng)
    r['Weighted mode'] = mode(bx, by, sex, sey, True, nboot, rng=rng)
    return r

def leave_one_out(bx, by, sey, ids):
    out = []
    for i in range(len(bx)):
        m = np.ones(len(bx), bool); m[i] = False
        r = ivw(bx[m], by[m], sey[m]); out.append((ids[i], r['b'], r['se'], r['pval']))
    return out

# ---------------- MR-PRESSO (Verbanck et al. 2018), following the R package ----------------
def _loo_rss(X, Y):
    # weighted data already; leave-one-out slope through origin; returns per-SNP squared residuals
    n = len(X); res = np.empty(n)
    sxx = np.sum(X**2); sxy = np.sum(X * Y)
    for i in range(n):
        b = (sxy - X[i] * Y[i]) / (sxx - X[i]**2)
        res[i] = (Y[i] - b * X[i])**2
    return res

def presso(bx, by, sex, sey, nsim=1000, signif=0.05, seed=1):
    rng = np.random.default_rng(seed)
    n = len(bx); wts = 1 / sey**2; sw = np.sqrt(wts)
    X = bx * sw; Y = by * sw
    rss_obs_i = _loo_rss(X, Y); rss_obs = rss_obs_i.sum()
    # simulate expected distribution under no horizontal pleiotropy
    sim_tot = np.empty(nsim); sim_i = np.empty((nsim, n))
    # leave-one-out slope for generating data
    b_loo = np.array([ (np.sum(X*Y) - X[i]*Y[i]) / (np.sum(X**2) - X[i]**2) for i in range(n)])  # weighted scale slope (same as unweighted slope)
    for s in range(nsim):
        bxs = rng.normal(bx, sex)
        bys = rng.normal(bx * b_loo, sey)
        Xs = bxs * sw; Ys = bys * sw
        r = _loo_rss(Xs, Ys); sim_i[s] = r; sim_tot[s] = r.sum()
    global_p = (np.sum(sim_tot > rss_obs) + 0) / nsim
    global_p = max(global_p, 1.0 / nsim)
    out_p = np.minimum(1.0, (np.sum(sim_i >= rss_obs_i[None, :], axis=0) / nsim) * n)
    outliers = np.where(out_p < signif)[0] if global_p < signif else np.array([], int)
    keep = np.setdiff1d(np.arange(n), outliers)
    raw = ivw(bx, by, sey)
    corr = ivw(bx[keep], by[keep], sey[keep]) if len(outliers) else raw
    return dict(global_rss=rss_obs, global_p=global_p, outlier_idx=outliers, outlier_p=out_p,
                raw=raw, corrected=corr, keep=keep)
