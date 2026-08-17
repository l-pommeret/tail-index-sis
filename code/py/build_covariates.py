"""Matrice de covariables, calculee sur LA MEME fenetre que la reponse.

Deux exigences que la premiere version violait.

1. Coherence X / Y. La reponse Y = exp(max NLL) porte sur une fenetre de
   NEED = 576 jetons. Calculer les covariables sur le document entier ferait
   porter X et Y sur deux objets differents. Tout est donc calcule sur la
   fenetre decodee, ce qui a un second benefice : les fenetres ont toutes la
   meme longueur en jetons (~2600 caracteres), donc la longueur cesse d'etre
   une nuisance d'echelle ET les frequences deviennent comparables d'un
   document a l'autre.

2. Densite des trigrammes. Le filtre de continuite ecarte toute colonne portant
   un atome de plus de MAX_ATOM = 2 %. Un trigramme present dans une fraction
   df des documents porte un atome de zero de 1 - df : il faut donc df >= 0.98
   pour passer, pas 0.30. Sur fenetres de 576 jetons, df >= 0.80 laisse ~200
   trigrammes et df >= 0.60 environ 450 ; le filtre tranche ensuite.

Les trigrammes sont choisis sur leur frequence documentaire seule, jamais sur
leur lien avec Y : un choix guide par la reponse invaliderait le criblage.

usage: python code/py/build_covariates.py DOCS.parquet OUT.npz [N] [NPROC] [DFMIN]
"""
from __future__ import annotations

import collections
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")
os.environ.setdefault("HF_HOME", "/people/pommeret/.cache/huggingface")
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from doc_covariates import surface_features, trigram_counter, trigram_features

DOCS = sys.argv[1] if len(sys.argv) > 1 else "results/wild/main_docs.parquet"
OUT = sys.argv[2] if len(sys.argv) > 2 else "results/wild/cov.npz"
NDOC = int(sys.argv[3]) if len(sys.argv) > 3 else 5000
NPROC = int(sys.argv[4]) if len(sys.argv) > 4 else 30
DF_MIN = float(sys.argv[5]) if len(sys.argv) > 5 else 0.60

MODEL = "EleutherAI/pythia-410m"
NEED = 576                 # BURN + T de extract_nll.py
MIN_DISTINCT = 100
MAX_ATOM = 0.02
NTRI_MAX = 1200


def _surface_chunk(texts):
    return [surface_features(t) for t in texts]


def _tri_chunk(args):
    texts, keys = args
    return np.stack([trigram_features(t, keys) for t in texts])


def main() -> None:
    import multiprocessing as mp
    from transformers import AutoTokenizer

    df = pd.read_parquet(DOCS, columns=["doc_id", "text"])
    tok = AutoTokenizer.from_pretrained(MODEL)

    # --- fenetrage : exactement les jetons que extract_nll.py voit ------------
    t0 = time.time()
    win, ids = [], []
    B = 2000
    for s in range(0, len(df), B):
        part = df.iloc[s:s + B]
        enc = tok(part["text"].tolist(), truncation=True,
                  max_length=NEED)["input_ids"]
        for e, did in zip(enc, part["doc_id"].values):
            if len(e) >= NEED:
                win.append(tok.decode(e)); ids.append(int(did))
        if len(win) >= NDOC:
            break
    win, ids = win[:NDOC], ids[:NDOC]
    n = len(win)
    L = np.array([len(w) for w in win], dtype=float)
    print(f"{n} fenetres de {NEED} jetons ({time.time()-t0:.0f} s)")
    print(f"  caracteres par fenetre : med {np.median(L):.0f}, "
          f"q10 {np.quantile(L,.1):.0f}, q90 {np.quantile(L,.9):.0f}")

    # --- choix des trigrammes, sans regarder Y -------------------------------
    t0 = time.time()
    dfreq = collections.Counter()
    for t in win:
        dfreq.update(set(trigram_counter(t)))
    keys = [k for k, v in dfreq.most_common() if v >= DF_MIN * n][:NTRI_MAX]
    print(f"  trigrammes candidats a df >= {DF_MIN:.2f} : {len(keys)} "
          f"({time.time()-t0:.0f} s)")

    chunks = [win[i:i + 250] for i in range(0, n, 250)]
    with mp.Pool(NPROC) as pool:
        t0 = time.time()
        surf = [r for part in pool.imap(_surface_chunk, chunks) for r in part]
        tri = np.concatenate(
            pool.map(_tri_chunk, [(c, keys) for c in chunks]), axis=0)
    print(f"  covariables calculees en {time.time()-t0:.0f} s")

    S = pd.DataFrame(surf).astype(np.float64)
    X = np.concatenate([S.values, tri.astype(np.float64)], axis=1)
    names = list(S.columns) + [f"tri_{i:04d}" for i in range(tri.shape[1])]
    tri_keys = [""] * len(S.columns) + keys
    block = ["surface"] * len(S.columns) + ["trigram"] * tri.shape[1]

    # --- filtre de continuite ------------------------------------------------
    keep = np.ones(X.shape[1], dtype=bool)
    for j in range(X.shape[1]):
        v, c = np.unique(X[:, j], return_counts=True)
        keep[j] = (len(v) >= MIN_DISTINCT) and (c.max() / n <= MAX_ATOM)
    ns = sum(1 for j in range(len(S.columns)) if keep[j])
    nt = keep.sum() - ns
    print(f"  filtre de continuite : {keep.sum()} retenues sur {len(keep)} "
          f"({ns} surface + {nt} trigrammes)")
    drop = [names[j] for j in range(len(S.columns)) if not keep[j]]
    if drop:
        print(f"    surface ecartees ({len(drop)}) : {', '.join(drop)}")

    X = X[:, keep]
    names = [names[j] for j, k in enumerate(keep) if k]
    tri_keys = [tri_keys[j] for j, k in enumerate(keep) if k]
    block = [block[j] for j, k in enumerate(keep) if k]
    os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
    np.savez_compressed(OUT, X=X.astype(np.float32),
                        names=np.array(names, dtype=object),
                        tri_keys=np.array(tri_keys, dtype=object),
                        block=np.array(block, dtype=object),
                        doc_id=np.array(ids, dtype=np.int64))
    # les fenetres servent aussi a l'encodeur : on les conserve
    pd.DataFrame({"doc_id": ids, "window": win}).to_parquet(
        OUT.replace(".npz", "_windows.parquet"), index=False)
    print(f"ECRIT {OUT} : {X.shape[0]} x {X.shape[1]}")


if __name__ == "__main__":
    main()
