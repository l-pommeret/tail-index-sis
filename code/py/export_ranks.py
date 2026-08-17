"""Pont binaire entre les activations (parquet) et le criblage (R).

Le score ne depend des covariables que par leurs RANGS, donc on transforme en
rangs cote Python et on ecrit une matrice binaire compacte plutot qu'un CSV de
plusieurs gigaoctets. Le departage des ex aequo est aleatoire a graine fixee,
comme dans screen_llm.R : la geometrie des fenetres en depend.

Configurations produites :
  cN     une seule couche N (p = 1024) -- la barre de competition est la plus
         basse et la question devient "ou dans la profondeur ?"
  c3     trois couches separees, 6/12/18 (p = 3072)
  acp    ACP sur les huit couches concatenees, K composantes -- meme nombre de
         coordonnees mais base tournee : repond a "existe-t-il UNE direction ?"
         quand le brut repond a "existe-t-il un NEURONE ?"

usage: python code/py/export_ranks.py PARQUET OUTDIR CONFIG [K]
"""
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

PARQ = sys.argv[1]
OUTDIR = sys.argv[2]
CONFIG = sys.argv[3]
K = int(sys.argv[4]) if len(sys.argv) > 4 else 1024
SEED = 20260815

os.makedirs(OUTDIR, exist_ok=True)
df = pd.read_parquet(PARQ)
y = df["max_S_suffix"].values.astype(np.float64)
n = len(df)
acts = [c for c in df.columns if c.startswith("h")]
layers = sorted({c[1:3] for c in acts})
print(f"{n} documents, {len(acts)} activations, couches {layers}")

if CONFIG == "c3":                      # teste AVANT le motif "cNN", sinon
    cols = [c for c in acts if c[1:3] in ("06", "12", "18")]   # "c3" serait lu
elif CONFIG.startswith("c") and CONFIG[1:].isdigit():          # comme la couche 3
    lay = f"{int(CONFIG[1:]):02d}"
    cols = [c for c in acts if c[1:3] == lay]
elif CONFIG == "acp":
    cols = acts
else:
    raise SystemExit(f"config inconnue : {CONFIG}")

X = df[cols].values.astype(np.float64)
print(f"config {CONFIG} : {X.shape[1]} colonnes brutes")

if CONFIG == "acp":
    Xc = X - X.mean(0, keepdims=True)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    k = min(K, Vt.shape[0])
    X = Xc @ Vt[:k].T
    var = (S ** 2) / (S ** 2).sum()
    print(f"ACP {k} composantes, variance expliquee {var[:k].sum():.4f}")
    cols = [f"pc_{i:04d}" for i in range(k)]

# rangs, ex aequo departages aleatoirement a graine fixee
rng = np.random.default_rng(SEED)
R = np.empty_like(X)
for j in range(X.shape[1]):
    noise = rng.random(n) * 1e-12 * (np.abs(X[:, j]).max() + 1e-12)
    R[:, j] = np.argsort(np.argsort(X[:, j] + noise)) + 1.0

base = os.path.join(OUTDIR, CONFIG)
np.ascontiguousarray(R.astype(np.float64).T).tofile(base + "_ranks.bin")
y.astype(np.float64).tofile(base + "_y.bin")
with open(base + "_meta.txt", "w") as f:
    f.write(f"{n}\n{R.shape[1]}\n")
    f.write("\n".join(cols) + "\n")
print(f"ECRIT {base}_ranks.bin ({R.nbytes/1e9:.2f} Go), {R.shape[1]} coordonnees")
print(f"seuil de competition a p={R.shape[1]} : "
      f"{abs(np.percentile(np.random.default_rng(0).normal(size=200000), 100/R.shape[1])):.2f} sigma")
