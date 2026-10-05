"""Minimal Python port of limma's two-group moderated t-test (lmFit + eBayes).

Follows Smyth (2004), Stat Appl Genet Mol Biol 3:3, and limma's fitFDist / trigammaInverse.
Only the two-group, no-covariate design is needed here.
"""
import numpy as np
import pandas as pd
from scipy import special, stats
from statsmodels.stats.multitest import multipletests


def trigamma_inverse(x):
    """Solve trigamma(y) = x for y (Newton iteration, as in limma::trigammaInverse)."""
    x = np.asarray(x, dtype=float)
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if np.all(-dif / y < 1e-8):
            break
    return y


def fit_f_dist(s2, df):
    """Estimate prior d0 and s0^2 from residual variances (limma::fitFDist, no covariate)."""
    s2 = np.maximum(s2, 1e-5 * np.median(s2))
    z = np.log(s2)
    e = z - special.digamma(df / 2) + np.log(df / 2)
    emean = e.mean()
    evar = ((e - emean) ** 2).sum() / (len(e) - 1) - special.polygamma(1, df / 2)
    if evar > 0:
        d0 = 2 * trigamma_inverse(evar)
        s02 = np.exp(emean + special.digamma(d0 / 2) - np.log(d0 / 2))
    else:
        d0, s02 = np.inf, np.exp(emean)
    return float(d0), float(s02)


def moderated_t(expr: pd.DataFrame, case: list, control: list) -> pd.DataFrame:
    """expr: genes x samples (log2). Returns log2FC (case - control), moderated t, p, BH adj p."""
    a = expr[case].to_numpy(float)
    b = expr[control].to_numpy(float)
    na, nb = a.shape[1], b.shape[1]
    lfc = a.mean(axis=1) - b.mean(axis=1)
    df = na + nb - 2
    s2 = (((a - a.mean(1, keepdims=True)) ** 2).sum(1) + ((b - b.mean(1, keepdims=True)) ** 2).sum(1)) / df
    d0, s02 = fit_f_dist(s2, df)
    if np.isinf(d0):
        s2_post, df_total = np.full_like(s2, s02), np.inf
    else:
        s2_post = (d0 * s02 + df * s2) / (d0 + df)
        df_total = min(df + d0, df * len(s2))  # limma caps at the pooled residual df
    se = np.sqrt(s2_post * (1 / na + 1 / nb))
    t = lfc / se
    p = 2 * (stats.norm.sf(np.abs(t)) if np.isinf(df_total) else stats.t.sf(np.abs(t), df_total))
    adj = multipletests(p, method="fdr_bh")[1]
    out = pd.DataFrame({"log2FC": lfc, "t": t, "p_value": p, "adj_p": adj,
                        "mean_expr": expr[case + control].mean(axis=1).to_numpy()}, index=expr.index)
    out.attrs.update(d0=d0, s02=s02, df_residual=df)
    return out
