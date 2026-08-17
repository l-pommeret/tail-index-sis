"""Assemblage X (surface + trigrammes + encodeur) et Y, avec export vers R.

Le score ne depend des covariables que par leurs rangs : on transforme donc en
rangs cote Python et on ecrit une matrice binaire, comme export_ranks.py. Le
departage des ex aequo est aleatoire a graine fixee -- la geometrie des
fenetres en depend, et un departage par indice de colonne ferait fuir l'ordre
des covariables dans le resultat (cf. l'audit RERUN_RANKFIX).

usage: python code/py/assemble_wild.py OUTDIR
"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "results/wild/ranks"
SEED = 20260817
os.makedirs(OUTDIR, exist_ok=True)

c = np.load("results/wild/cov.npz", allow_pickle=True)
e = np.load("results/wild/enc.npz", allow_pickle=True)
y = pd.read_parquet("results/wild/main_nll.parquet")

# alignement strict par doc_id : les trois etages ont filtre independamment
ids = np.intersect1d(np.intersect1d(c["doc_id"], e["doc_id"]), y["doc_id"].values)
ic = np.searchsorted(c["doc_id"], ids, sorter=np.argsort(c["doc_id"]))
ic = np.argsort(c["doc_id"])[ic]
ie = np.argsort(e["doc_id"])[np.searchsorted(e["doc_id"], ids,
                                             sorter=np.argsort(e["doc_id"]))]
yy = y.set_index("doc_id").loc[ids]
assert (c["doc_id"][ic] == ids).all() and (e["doc_id"][ie] == ids).all()

X = np.concatenate([c["X"][ic].astype(np.float64),
                    e["E"][ie].astype(np.float64)], axis=1)
names = list(c["names"]) + list(e["names"])
block = list(c["block"]) + ["encoder"] * e["E"].shape[1]
n, p = X.shape
print(f"n = {n}, p = {p} ({block.count('surface')} surface + "
      f"{block.count('trigram')} trigrammes + {block.count('encoder')} encodeur)")
print(f"  n*alpha a alpha=0.1 : {0.1*n:.0f}  -> p/n*alpha = {p/(0.1*n):.2f}")

M = yy["max_nll"].values.astype(np.float64)     # log Y
rng = np.random.default_rng(SEED)
R = np.empty_like(X)
for j in range(p):
    noise = rng.random(n) * 1e-12 * (np.abs(X[:, j]).max() + 1e-12)
    R[:, j] = np.argsort(np.argsort(X[:, j] + noise)) + 1.0

R.tofile(os.path.join(OUTDIR, "wild_ranks.bin"))
M.tofile(os.path.join(OUTDIR, "wild_y.bin"))
with open(os.path.join(OUTDIR, "wild_meta.txt"), "w") as f:
    f.write(f"{n}\n{p}\n" + "\n".join(names) + "\n")
pd.DataFrame({"name": names, "block": block}).to_csv(
    os.path.join(OUTDIR, "wild_blocks.csv"), index=False)
print(f"ECRIT {OUTDIR}/wild_ranks.bin ({R.nbytes/1e6:.0f} Mo)")
print(f"  barre de competition a p={p} : "
      f"{abs(np.percentile(rng.normal(size=400000), 100/p)):.2f} sigma")
