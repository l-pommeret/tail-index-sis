"""Coordonnees d'encodeur : le fond de nuisance dense.

Role dans le design. Les covariables de surface interpretables sont peu
nombreuses une fois le filtre de continuite passe (~35 sur 63 : les taux de
ponctuation rare sont gonfles de zeros sur une fenetre de 2600 caracteres).
Les coordonnees d'encodeur sont continues par construction et jouent le fond
de nuisance correle des simulations -- l'analogue de l'AR(1) de la famille A.

Elles sont calculees sur LA MEME fenetre que la reponse et que la surface.

usage: python code/py/encode_docs.py WINDOWS.parquet OUT.npz [MODEL] [GPU]
"""
import os, sys, time, warnings; warnings.filterwarnings("ignore")
os.environ.setdefault("HF_HOME", "/people/pommeret/.cache/huggingface")
GPU = sys.argv[4] if len(sys.argv) > 4 else "1"
os.environ["CUDA_VISIBLE_DEVICES"] = GPU
import numpy as np, pandas as pd

WIN = sys.argv[1] if len(sys.argv) > 1 else "results/wild/cov_windows.parquet"
OUT = sys.argv[2] if len(sys.argv) > 2 else "results/wild/enc.npz"
MODEL = sys.argv[3] if len(sys.argv) > 3 else "sentence-transformers/all-MiniLM-L6-v2"

from sentence_transformers import SentenceTransformer
d = pd.read_parquet(WIN)
m = SentenceTransformer(MODEL, device="cuda")
t0 = time.time()
E = m.encode(d["window"].tolist(), batch_size=128, convert_to_numpy=True,
             normalize_embeddings=False, show_progress_bar=False)
print(f"{MODEL} : {E.shape[0]} x {E.shape[1]} en {time.time()-t0:.0f} s")

# controle de continuite : une coordonnee d'encodeur doit etre continue
n = E.shape[0]
atom = np.array([np.unique(E[:, j], return_counts=True)[1].max() / n
                 for j in range(E.shape[1])])
ndist = np.array([len(np.unique(E[:, j])) for j in range(E.shape[1])])
keep = (ndist >= 100) & (atom <= 0.02)
print(f"  filtre de continuite : {keep.sum()} sur {E.shape[1]} "
      f"(atome max {atom.max():.4f})")
np.savez_compressed(OUT, E=E[:, keep].astype(np.float32),
                    names=np.array([f"enc_{j:04d}" for j in np.where(keep)[0]],
                                   dtype=object),
                    doc_id=d["doc_id"].values)
print(f"ECRIT {OUT} : {E[:, keep].shape}")
