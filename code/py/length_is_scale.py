"""La longueur agit-elle sur l'echelle et jamais sur l'indice ?

Pour T variables a variation reguliere d'indice -1/gamma, faiblement
dependantes d'indice extremal theta,

    P(max_{t<=T} > y | x) ~ theta T c(x) y^{-1/gamma(x)},

donc les quantiles de Y sont multiplies par T^gamma et l'indice gamma est
inchange. En passant au logarithme, avec M = max_t S_t = log Y,

    E[M | x] = gamma(x) log T + O(1),      Var[M | x] = pi^2 gamma(x)^2 / 6.

Trois consequences testables, et c'est le point : le statut de la longueur est
etabli analytiquement, pas mesure. Sur des donnees reelles on peut donc
verifier la theorie plutot que la supposer.

  (a) la mediane de M croit affinement en log T, de pente gamma ;
  (b) gamma estime par strate de longueur est constant ;
  (c) l'ecart-type de M par strate vaut environ pi gamma / sqrt(6) = 1.283 gamma.

usage: python code/py/length_is_scale.py LIBRE.parquet [FIXE.parquet]
"""
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

sys.path.insert(0, "code/py")
from diagnostics import fit_gpd, hill, profile_ci_xi

LIBRE = sys.argv[1] if len(sys.argv) > 1 else "results/llm/pilot_libre.parquet"
FIXE = sys.argv[2] if len(sys.argv) > 2 else "results/llm/pilot_fixe.parquet"
NQ = 5


def strata_report(df, by, label, out):
    qs = np.quantile(df[by], np.linspace(0, 1, NQ + 1))
    out.append(f"\n  strates par {label} :")
    out.append(f"    {'strate':>7s} {'n':>5s} {'T med':>8s} {'log T':>7s} "
               f"{'med M':>7s} {'sd M':>6s} {'Hill':>7s} {'xi':>7s} {'IC95':>18s}")
    rows = []
    for i in range(NQ):
        lo, hi = qs[i], qs[i + 1]
        sel = (df[by] >= lo) & (df[by] <= hi) if i == NQ - 1 else \
              (df[by] >= lo) & (df[by] < hi)
        d = df[sel]
        if len(d) < 120:
            continue
        M = np.sort(d.max_S_small.values)[::-1]
        k = int(0.25 * len(M))
        u = M[k]
        exc = np.exp(M[:k] - u) - 1.0
        xi, ci = profile_ci_xi(exc)
        h = hill(M, k)
        Tmed = float(d.T_scored.median())
        rows.append((np.log(Tmed), float(np.median(d.max_S_small)), h, xi))
        out.append(f"    {'Q'+str(i+1):>7s} {len(d):5d} {Tmed:8.0f} "
                   f"{np.log(Tmed):7.3f} {np.median(d.max_S_small):7.3f} "
                   f"{np.std(d.max_S_small):6.3f} {h:7.4f} {xi:+7.4f} "
                   f"[{ci[0]:+.3f},{ci[1]:+.3f}]")
    return rows


def main():
    out = ["LA LONGUEUR EST-ELLE UNE NUISANCE D'ECHELLE ?", ""]
    lib = pd.read_parquet(LIBRE)
    fix = pd.read_parquet(FIXE)

    out.append(f"mode libre : n={len(lib)}, T de {lib.T_scored.min():.0f} a "
               f"{lib.T_scored.max():.0f} (med {lib.T_scored.median():.0f})")
    out.append(f"mode fixe  : n={len(fix)}, T = {fix.T_scored.iloc[0]:.0f} constant")

    out.append("\n=== (a) et (b) : mediane de M et indice, par strate de longueur ===")
    rows = strata_report(lib, "T_scored", "T_scored (mode libre)", out)

    if len(rows) >= 3:
        lt = np.array([r[0] for r in rows])
        med = np.array([r[1] for r in rows])
        hills = np.array([r[2] for r in rows])
        xis = np.array([r[3] for r in rows])
        A = np.vstack([lt, np.ones_like(lt)]).T
        slope, icept = np.linalg.lstsq(A, med, rcond=None)[0]
        resid = med - (A @ [slope, icept])
        out.append(f"\n  (a) regression de med(M) sur log T : pente = {slope:.3f} nats/nat, "
                   f"ordonnee {icept:.3f}, residu max {np.abs(resid).max():.3f}")
        out.append(f"      gamma estime par ailleurs : Hill moyen {hills.mean():.3f}, "
                   f"xi moyen {xis.mean():.3f}")
        out.append(f"      la theorie predit pente = gamma. Ecart pente - Hill = "
                   f"{slope - hills.mean():+.3f}")
        out.append(f"\n  (b) gamma par strate : Hill de {hills.min():.3f} a {hills.max():.3f} "
                   f"(etendue {np.ptp(hills):.3f}), xi de {xis.min():.3f} a {xis.max():.3f}")
        out.append(f"      a comparer a la variation de med(M) sur les memes strates : "
                   f"{np.ptp(med):.3f} nats")

    out.append("\n=== (c) dispersion de M par strate ===")
    sds = []
    qs = np.quantile(lib.T_scored, np.linspace(0, 1, NQ + 1))
    for i in range(NQ):
        lo, hi = qs[i], qs[i + 1]
        sel = (lib.T_scored >= lo) & (lib.T_scored <= hi) if i == NQ - 1 else \
              (lib.T_scored >= lo) & (lib.T_scored < hi)
        if sel.sum() >= 120:
            sds.append(float(np.std(lib.max_S_small[sel])))
    if sds:
        g_from_sd = np.mean(sds) * np.sqrt(6) / np.pi
        out.append(f"  sd(M) par strate : {', '.join(f'{s:.3f}' for s in sds)}")
        out.append(f"  gamma deduit de sd = pi gamma / sqrt(6) : {g_from_sd:.3f}")
        out.append("  (borne basse attendue : la variance inter-documents s'ajoute a"
                   " celle de Gumbel)")

    out.append("\n=== comparaison des deux modes ===")
    for nm, d in (("libre", lib), ("fixe", fix)):
        M = np.sort(d.max_S_small.values)[::-1]
        k = int(0.20 * len(M))
        out.append(f"  {nm:6s} : med(M) {np.median(d.max_S_small):6.3f}, "
                   f"Hill(20%) {hill(M, k):.4f}, "
                   f"rho(M, len_words) {pd.Series(d.len_words).corr(pd.Series(d.max_S_small), method='spearman'):+.3f}")
    out.append("\n  lecture : l'indice doit etre le meme dans les deux modes,"
               "\n  la correlation avec la longueur ne doit survivre qu'en mode libre.")

    txt = "\n".join(out)
    print(txt)
    with open("results/llm/length_is_scale.txt", "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
