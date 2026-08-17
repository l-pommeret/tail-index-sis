"""Export des covariables de perplexite vers le criblage R, en rangs.

Le filtre de continuite est reapplique ici : une NLL est continue, mais un
quantile bas d'un ecart entre modeles peut saturer. On verifie plutot que
supposer.

usage: python code/py/export_ppl_ranks.py IN.parquet OUTDIR
"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

IN = sys.argv[1] if len(sys.argv) > 1 else "results/wild/ppl/ppl_cov.parquet"
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else "results/wild/ppl_ranks"
SEED = 20260817
os.makedirs(OUTDIR, exist_ok=True)

d = pd.read_parquet(IN)
meta = ["doc_id", "Ylog", "suffix_mean"]
cols = [c for c in d.columns if c not in meta]
X = d[cols].values.astype(np.float64)
n, p0 = X.shape

keep = np.ones(p0, dtype=bool)
for j in range(p0):
    v, c = np.unique(X[:, j], return_counts=True)
    keep[j] = (len(v) >= 100) and (c.max() / n <= 0.02)
print(f"filtre de continuite : {keep.sum()} sur {p0} retenues")
if (~keep).any():
    print("  ecartees :", ", ".join(np.array(cols)[~keep][:10]))
X = X[:, keep]; cols = list(np.array(cols)[keep]); p = X.shape[1]

block = [c.split("_")[0] for c in cols]
print(f"n = {n}, p = {p}")
for b in sorted(set(block)):
    print(f"    {b:14s} {block.count(b)}")

M = d["Ylog"].values.astype(np.float64)
rng = np.random.default_rng(SEED)
R = np.empty_like(X)
for j in range(p):
    noise = rng.random(n) * 1e-12 * (np.abs(X[:, j]).max() + 1e-12)
    R[:, j] = np.argsort(np.argsort(X[:, j] + noise)) + 1.0

R.tofile(os.path.join(OUTDIR, "wild_ranks.bin"))
M.tofile(os.path.join(OUTDIR, "wild_y.bin"))
with open(os.path.join(OUTDIR, "wild_meta.txt"), "w") as f:
    f.write(f"{n}\n{p}\n" + "\n".join(cols) + "\n")
pd.DataFrame({"name": cols, "block": block}).to_csv(
    os.path.join(OUTDIR, "wild_blocks.csv"), index=False)
print(f"ECRIT {OUTDIR} ({R.nbytes/1e6:.0f} Mo)")
