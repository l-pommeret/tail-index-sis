"""Y = exp(max NLL) est-il a variation reguliere ? C'est le verrou du projet.

Si la surprisal a une queue superieure exponentielle de taux lambda, alors
Y = e^S est Pareto d'indice gamma = 1/lambda. Comme Hill sur Y vaut

    gamma_k = (1/k) sum_{i<=k} [ log Y_(i) - log Y_(k+1) ]

et que log Y = max_nll, l'estimateur est exactement l'exces moyen de max_nll
au-dessus de son (k+1)-eme plus grand : on travaille directement sur max_nll.

Trois diagnostics :
  1. plateau de Hill    gamma stable sur une plage de k = signature Pareto
  2. QQ exponentiel     log Y contre les quantiles exponentiels : droite
  3. IC de profil       via l'ecart-type asymptotique gamma/sqrt(k)

usage: python code/py/evt_gate.py NLL.parquet [COL]
"""
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

P = sys.argv[1] if len(sys.argv) > 1 else "results/wild/pilot_nll.parquet"
COL = sys.argv[2] if len(sys.argv) > 2 else "max_nll"

d = pd.read_parquet(P)
L = np.sort(d[COL].values)[::-1]          # log Y decroissant
n = len(L)
print(f"{P} : n = {n}, reponse log Y = {COL}")
print(f"  log Y : med {np.median(L):.3f}, q95 {np.quantile(L,.95):.3f}, "
      f"max {L.max():.3f}\n")

print("plateau de Hill (gamma = exces moyen de log Y au-dessus du seuil)")
print("       k   k/n    gamma    e.t.        IC95            seuil")
rows = []
for frac in (0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50):
    k = int(frac * n)
    if k < 20:
        continue
    g = float(np.mean(L[:k] - L[k]))
    se = g / np.sqrt(k)
    rows.append((k, frac, g, se))
    print(f"  {k:6d}  {frac:.2f}  {g:7.4f}  {se:.4f}  "
          f"[{g-1.96*se:.4f},{g+1.96*se:.4f}]  {L[k]:.3f}")

g = np.array([r[2] for r in rows])
print(f"\n  etendue du plateau sur k/n in [0.02,0.50] : "
      f"{g.min():.4f} a {g.max():.4f} (rapport {g.max()/g.min():.2f})")
sub = [r[2] for r in rows if 0.05 <= r[1] <= 0.25]
print(f"  gamma median sur k/n in [0.05,0.25] : {np.median(sub):.4f}")

# --- QQ exponentiel sur les exces -------------------------------------------
print("\nQQ exponentiel des exces (droite <=> queue exponentielle de log Y)")
for frac in (0.05, 0.10, 0.20):
    k = int(frac * n)
    e = L[:k] - L[k]
    # quantiles exponentiels croissants, alignes sur sort(e) croissant
    q = -np.log(1 - (np.arange(1, k + 1) - 0.5) / k)
    r = np.corrcoef(q, np.sort(e))[0, 1]
    # ecart relatif max a la droite des moindres carres
    b = np.polyfit(q, np.sort(e), 1)
    res = np.abs(np.sort(e) - np.polyval(b, q)).max() / e.max()
    print(f"  k/n = {frac:.2f} (k={k:5d}) : r = {r:.5f}, "
          f"ecart max relatif {res:.3f}, pente {b[0]:.4f}")

# --- comparaison a une alternative a queue fine ------------------------------
print("\ncontrole : meme diagnostic sur mean_nll (ne doit PAS etre Pareto)")
if "mean_nll" in d.columns:
    M = np.sort(d["mean_nll"].values)[::-1]
    for frac in (0.05, 0.20):
        k = int(frac * n)
        print(f"  k/n = {frac:.2f} : gamma = {np.mean(M[:k]-M[k]):.4f}")
