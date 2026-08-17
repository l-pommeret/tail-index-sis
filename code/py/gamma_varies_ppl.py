"""gamma(x) varie-t-il ? Diagnostic independant de l'ecran.

Un resultat de criblage nul a deux causes possibles : l'ecran manque de
puissance, ou gamma(x) ne varie pas et il n'y a rien a trouver. On tranche
sans passer par l'ecran.

On decoupe en Q strates de tailles egales selon chaque covariable et on estime
gamma dans chacune par Hill. La dispersion inter-strates melange la variation
reelle de gamma et le bruit de Hill (~ gamma/sqrt(k)). On calibre donc par des
strates ALEATOIRES de memes tailles, qui ne contiennent que le bruit.

usage: python code/py/gamma_varies_ppl.py PARQUET [B] [Q] [FRAC]
"""
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

P = sys.argv[1] if len(sys.argv) > 1 else "results/wild/ppl_same/ppl_cov.parquet"
B = int(sys.argv[2]) if len(sys.argv) > 2 else 500
Q = int(sys.argv[3]) if len(sys.argv) > 3 else 4
FRAC = float(sys.argv[4]) if len(sys.argv) > 4 else 0.20

d = pd.read_parquet(P)
M = d["Ylog"].values.astype(float)      # log Y : Hill = exces moyen
n = len(M)
cands = [c for c in d.columns if c not in ("doc_id", "Ylog", "suffix_mean")
         and pd.api.types.is_numeric_dtype(d[c]) and d[c].nunique() >= 100]
print(f"{P}\nn = {n}, {len(cands)} covariables, Q = {Q} strates, k/strate "
      f"= {int(FRAC*n/Q)}")

def sd_inter(groups):
    g = []
    for i in range(Q):
        m = np.sort(M[groups == i])[::-1]
        k = int(FRAC * len(m))
        if k < 30:
            return None
        g.append(np.mean(m[:k] - m[k]))
    return float(np.std(g, ddof=1)), g

rng = np.random.default_rng(7)
null = []
for _ in range(B):
    gr = rng.permutation(np.repeat(np.arange(Q), int(np.ceil(n/Q)))[:n])
    r = sd_inter(gr)
    if r: null.append(r[0])
null = np.array(null)
print(f"\nsous H0 ({len(null)} decoupages aleatoires) : med {np.median(null):.4f}, "
      f"q95 {np.quantile(null,.95):.4f}, q99 {np.quantile(null,.99):.4f}")
gall = np.mean(np.sort(M)[::-1][:int(FRAC*n)] - np.sort(M)[::-1][int(FRAC*n)])
print(f"  gamma global = {gall:.4f} ; prediction gamma/sqrt(k) = "
      f"{gall/np.sqrt(int(FRAC*n/Q)):.4f}")

rows = []
for v in cands:
    x = d[v].values
    qs = np.quantile(x, np.linspace(0, 1, Q+1))
    gr = np.clip(np.searchsorted(qs[1:-1], x, side="right"), 0, Q-1)
    if np.bincount(gr, minlength=Q).min() < n/(3*Q):
        continue
    r = sd_inter(gr)
    if not r: continue
    s, g = r
    rows.append((v, s, float(np.mean(null >= s)), min(g), max(g)))

res = pd.DataFrame(rows, columns=["covariate","sd_inter","p_perm","gmin","gmax"])
res = res.sort_values("sd_inter", ascending=False)
print(f"\n{len(res)} covariables evaluees")
print(f"  dispersion observee : med {res.sd_inter.median():.4f}, "
      f"max {res.sd_inter.max():.4f}  (H0 med {np.median(null):.4f})")
print(f"  p < 0.05 : {int((res.p_perm<.05).sum())} sur {len(res)} "
      f"(attendu {0.05*len(res):.0f})")
print(f"  p < 0.01 : {int((res.p_perm<.01).sum())} sur {len(res)} "
      f"(attendu {0.01*len(res):.0f})")
print("\n  dix plus fortes dispersions :")
print(res.head(10).to_string(index=False, float_format=lambda v: f"{v:.4f}"))
res.to_csv(P.replace(".parquet", "_gammavar.csv"), index=False)
