"""gamma(x) varie-t-il, tout court ?

Trois resultats nuls de criblage successifs peuvent avoir deux causes : soit
l'ecran manque de puissance, soit gamma(x) ne varie pas assez pour qu'il y ait
quelque chose a trouver. Cette question se tranche sans passer par l'ecran.

Protocole. Pour une covariable v, on decoupe en Q strates de tailles egales et
on estime gamma dans chacune par Hill. La dispersion inter-strates melange deux
choses : la variation reelle de gamma(x) et le bruit d'echantillonnage de Hill,
dont l'ecart-type vaut environ gamma / sqrt(k) avec k le nombre d'exces par
strate. On calibre donc par des strates ALEATOIRES de memes tailles : elles ne
contiennent que le bruit. Si le decoupage par covariable ne disperse pas plus
que le decoupage aleatoire, gamma ne depend pas de cette covariable.

usage: python code/py/gamma_varies.py PARQUET COLONNE_Y [B] [Q]
"""
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

sys.path.insert(0, "code/py")
from diagnostics import hill

PARQ = sys.argv[1] if len(sys.argv) > 1 else "results/llm/main_fixe.parquet"
YCOL = sys.argv[2] if len(sys.argv) > 2 else "max_S_small"
B = int(sys.argv[3]) if len(sys.argv) > 3 else 400
Q = int(sys.argv[4]) if len(sys.argv) > 4 else 4
FRAC = 0.20


def gamma_by_groups(M, groups, q):
    """Hill dans chaque groupe ; renvoie l'ecart-type inter-groupes."""
    g = []
    for i in range(q):
        m = np.sort(M[groups == i])[::-1]
        k = int(FRAC * len(m))
        if k < 30:
            return np.nan
        g.append(hill(m, k))
    return float(np.std(g, ddof=1)), g


def main():
    df = pd.read_parquet(PARQ)
    M = df[YCOL].values.astype(float)
    n = len(M)
    meta = {"doc_id", "url", "dump", "date", "argmax_tok", "top_S", "n_tok_doc",
            "T_scored", "max_S_small", "max_S_large", "max_delta",
            "mean_S_small", "mean_S_large", "max_S_suffix"}
    cands = [c for c in df.columns if c not in meta and
             pd.api.types.is_numeric_dtype(df[c]) and df[c].nunique() >= 100]
    print(f"{PARQ}\nn = {n}, reponse = {YCOL}, {len(cands)} covariables candidates")

    # --- loi nulle : strates aleatoires de memes tailles ---------------------
    rng = np.random.default_rng(7)
    null = []
    for _ in range(B):
        gr = rng.permutation(np.repeat(np.arange(Q), int(np.ceil(n / Q)))[:n])
        r = gamma_by_groups(M, gr, Q)
        if not isinstance(r, float):
            null.append(r[0])
    null = np.array(null)
    k_par_strate = int(FRAC * n / Q)
    print(f"\nloi nulle sur {len(null)} decoupages aleatoires en {Q} strates "
          f"({k_par_strate} exces par strate)")
    print(f"  ecart-type inter-strates sous H0 : med {np.median(null):.4f}, "
          f"q95 {np.quantile(null, .95):.4f}, q99 {np.quantile(null, .99):.4f}")
    print(f"  prediction theorique gamma/sqrt(k) : "
          f"{np.mean([hill(np.sort(M)[::-1], int(FRAC*n))]) / np.sqrt(k_par_strate):.4f}")

    # --- observe : strates par covariable ------------------------------------
    rows = []
    for v in cands:
        x = df[v].values
        qs = np.quantile(x, np.linspace(0, 1, Q + 1))
        gr = np.clip(np.searchsorted(qs[1:-1], x, side="right"), 0, Q - 1)
        if min(np.bincount(gr, minlength=Q)) < n / (3 * Q):
            continue
        r = gamma_by_groups(M, gr, Q)
        if isinstance(r, float):
            continue
        sd, g = r
        rows.append((v, sd, float(np.mean(null >= sd)), min(g), max(g)))

    res = pd.DataFrame(rows, columns=["covariate", "sd_inter", "p_perm",
                                      "gamma_min", "gamma_max"])
    res = res.sort_values("sd_inter", ascending=False)
    print(f"\n{len(res)} covariables evaluees ; dispersion inter-strates de gamma")
    print(f"  observee : med {res.sd_inter.median():.4f}, max {res.sd_inter.max():.4f}")
    print(f"  attendue sous H0 : med {np.median(null):.4f}")
    print(f"  p < 0.05 : {int((res.p_perm < .05).sum())} sur {len(res)} "
          f"(attendu par hasard {0.05*len(res):.0f})")
    print(f"  p < 0.01 : {int((res.p_perm < .01).sum())} sur {len(res)} "
          f"(attendu {0.01*len(res):.0f})")
    print("\n  dix covariables de plus forte dispersion :")
    print(res.head(10).to_string(index=False,
          float_format=lambda v: f"{v:.4f}"))
    res.to_csv(PARQ.replace(".parquet", f"_gammavar_{YCOL}.csv"), index=False)


if __name__ == "__main__":
    main()
