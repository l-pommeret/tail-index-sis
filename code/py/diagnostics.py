"""Diagnostics du pilote, avant tout screening.

Trois questions, dans l'ordre ou une reponse negative arrete tout :

1. PLAFOND (point 3). Y = exp(max_t S_t), donc log Y = max_t S_t. Si les plus
   grandes surprises s'accumulent contre une borne -- numerique ou d'ecart de
   logits -- la queue est tronquee et gamma est biaise vers le bas. On regarde
   la distribution des maxima et la part de documents dont le maximum tombe
   dans un intervalle etroit sous le maximum global.

2. GAMMA > 0 (points 1 et 6). C'est le risque central : si chaque loi
   conditionnelle est lognormale, gamma(x) = 0 partout, (C1) est violee et
   l'ecran classe du bruit. Le Hill plot ne suffit pas -- il ne teste rien. On
   ajuste une GPD par maximum de vraisemblance au-dessus de seuils successifs
   et on construit un IC de vraisemblance profilee sur xi. Comme log Y est un
   maximum de T surprises, l'estimateur de Hill de Y est exactement la fonction
   d'exces moyen de max_t S_t, ce qui rend la lecture directe.

3. GAMMA > 0 PAR STRATES. C'est la seule version conditionnelle du diagnostic,
   donc la seule qui teste vraiment (C1) : gamma marginal > 0 est compatible
   avec gamma(x) = 0 partout, par melange sur des types de documents.

usage: python code/py/diagnostics.py PILOT.parquet
"""
import json
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import genpareto

PATH = sys.argv[1] if len(sys.argv) > 1 else "results/llm/pilot_fixe.parquet"


# --------------------------------------------------------------- GPD par MLE --
def gpd_nll(params, x):
    xi, sigma = params
    if sigma <= 0:
        return np.inf
    z = x / sigma
    if xi >= 0:
        if np.any(1 + xi * z <= 0):
            return np.inf
        return len(x) * np.log(sigma) + (1 + 1 / xi) * np.sum(np.log1p(xi * z)) \
            if abs(xi) > 1e-8 else len(x) * np.log(sigma) + np.sum(z)
    if np.any(1 + xi * z <= 0):
        return np.inf
    return len(x) * np.log(sigma) + (1 + 1 / xi) * np.sum(np.log1p(xi * z))


def fit_gpd(x):
    best, bp = np.inf, None
    for xi0 in (-0.2, 0.05, 0.3, 0.8):
        for s0 in (np.mean(x), np.std(x) + 1e-9):
            r = minimize(gpd_nll, [xi0, max(s0, 1e-9)], args=(x,),
                         method="Nelder-Mead",
                         options=dict(maxiter=4000, xatol=1e-8, fatol=1e-10))
            if r.fun < best:
                best, bp = r.fun, r.x
    return bp[0], bp[1], best


def profile_ci_xi(x, level=0.95, grid=None):
    """IC de vraisemblance profilee sur xi : {xi : 2(l_max - l_prof(xi)) < q}."""
    from scipy.stats import chi2
    xi_hat, _, nll_min = fit_gpd(x)
    q = chi2.ppf(level, 1)
    if grid is None:
        grid = np.linspace(max(xi_hat - 0.6, -0.9), xi_hat + 0.9, 220)
    keep = []
    for xi in grid:
        r = minimize_scalar(lambda s: gpd_nll([xi, s], x),
                            bounds=(1e-9, 10 * (np.mean(x) + 1e-9)),
                            method="bounded",
                            options=dict(xatol=1e-10))
        if 2 * (r.fun - nll_min) < q:
            keep.append(xi)
    return xi_hat, (min(keep), max(keep)) if keep else (np.nan, np.nan)


def hill(m_sorted_desc, k):
    """Hill de Y a partir de M = log Y, deja trie par ordre decroissant.

    L'estimateur de Hill est la moyenne des log Y_(i) - log Y_(k+1) sur les k
    plus grandes, c'est-a-dire exactement la moyenne des M_(i) - M_(k+1) : ne
    pas reprendre de logarithme ici.
    """
    return float(np.mean(m_sorted_desc[:k] - m_sorted_desc[k]))


def report_tail(name, M, out):
    """M = log Y = max_t S_t. Hill sur Y = exces moyen de M."""
    Ms = np.sort(M)[::-1]
    n = len(Ms)
    out.append(f"\n--- {name} (n={n}) ---")
    out.append("  Hill (= exces moyen de max_t S) :")
    for frac in (0.02, 0.05, 0.10, 0.20):
        k = max(int(frac * n), 5)
        out.append(f"    k={k:5d} ({frac:.0%})  gamma_hat = {hill(Ms, k):.4f}")
    for frac in (0.10, 0.20):
        k = int(frac * n)
        u = Ms[k]
        exc = np.exp(Ms[:k] - u) - 1.0          # exces relatifs de Y au seuil e^u
        xi, ci = profile_ci_xi(exc)
        verdict = "gamma > 0" if ci[0] > 0 else ("indetermine" if ci[1] > 0 else "gamma = 0")
        out.append(f"  GPD au seuil {frac:.0%} (u={u:.2f} nats, {k} exces) : "
                   f"xi = {xi:.4f}, IC95 profile [{ci[0]:.4f}, {ci[1]:.4f}]  -> {verdict}")


def main():
    df = pd.read_parquet(PATH)
    out = [f"DIAGNOSTICS  {PATH}  n = {len(df)} documents",
           f"mode : T_scored median = {df.T_scored.median():.0f} tokens notes"]

    # ------------------------------------------------------- 1. plafond ------
    M = df.max_S_small.values
    top = np.sort(M)[::-1]
    out.append("\n=== 1. PLAFOND ===")
    out.append(f"max_t S : min {M.min():.2f}, med {np.median(M):.2f}, "
               f"q99 {np.quantile(M, .99):.2f}, max {M.max():.2f} nats")
    for w in (0.10, 0.25, 0.50, 1.00):
        share = float(np.mean(M > M.max() - w))
        out.append(f"  part des documents a moins de {w:.2f} nats du maximum global : {share:.4f}")
    gaps = np.diff(top[:20])
    out.append(f"  ecarts entre les 20 plus grands maxima : "
               f"med {np.median(-gaps):.3f}, min {(-gaps).min():.3f} nats")
    out.append("  lecture : une accumulation contre une borne se verrait comme une part"
               "\n  elevee juste sous le maximum et des ecarts qui s'effondrent.")

    # ------------------------------------------- 2. gamma > 0, marginal ------
    out.append("\n=== 2. GAMMA > 0 (marginal) ===")
    report_tail("Y = exp(max S) du petit modele", df.max_S_small.values, out)
    report_tail("Y = exp(max S) du grand modele", df.max_S_large.values, out)
    report_tail("Y_delta = exp(max (S_petit - S_grand))", df.max_delta.values, out)

    # ------------------------------------------- 3. gamma > 0 par strates ----
    out.append("\n=== 3. GAMMA > 0 PAR STRATES (test conditionnel de (C1)) ===")
    strat_vars = ["gzip_ratio", "rate_digit", "ttr", "byte_entropy"]
    for v in strat_vars:
        if v not in df.columns:
            continue
        qs = np.quantile(df[v], [0, .25, .5, .75, 1.0])
        out.append(f"\n  strates par {v} :")
        for i in range(4):
            lo, hi = qs[i], qs[i + 1]
            sel = (df[v] >= lo) & (df[v] <= hi) if i == 3 else (df[v] >= lo) & (df[v] < hi)
            Ms = np.sort(df.max_S_small.values[sel.values])[::-1]
            if len(Ms) < 80:
                continue
            k = int(0.20 * len(Ms))
            u = Ms[k]
            exc = np.exp(Ms[:k] - u) - 1.0
            xi, ci = profile_ci_xi(exc)
            out.append(f"    Q{i+1} (n={sel.sum():4d}) : Hill {hill(Ms, k):.4f}, "
                       f"xi {xi:+.4f}, IC95 [{ci[0]:+.4f}, {ci[1]:+.4f}]")

    # ------------------------- 4. la longueur agit-elle sur l'echelle seule ? --
    out.append("\n=== 4. LONGUEUR : echelle ou indice ? ===")
    if df.T_scored.nunique() > 1:
        out.append("  (mode libre) correlation de Spearman avec max_t S :")
        for v in ["len_chars", "len_words", "len_gzip", "T_scored"]:
            if v in df.columns:
                r = pd.Series(df[v]).corr(pd.Series(M), method="spearman")
                out.append(f"    {v:12s} rho = {r:+.3f}")
    else:
        out.append(f"  (mode fixe) T_scored constant = {int(df.T_scored.iloc[0])}, "
                   "aucun confondant de longueur par construction.")
        out.append("  correlations residuelles avec les proxys de longueur du document :")
        for v in ["len_chars", "len_words", "len_gzip"]:
            if v in df.columns:
                r = pd.Series(df[v]).corr(pd.Series(M), method="spearman")
                out.append(f"    {v:12s} rho = {r:+.3f}")

    txt = "\n".join(out)
    print(txt)
    with open(PATH.replace(".parquet", "_diagnostics.txt"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
