"""Covariables DERIVEES DE LA PERPLEXITE, et reponse sur un span disjoint.

Le probleme que ceci resout. Les covariables de surface (taux de ponctuation,
trigrammes) sont gonflees de zeros sur une fenetre de 2600 caracteres : 28 des
63 covariables de surface et 411 des 449 trigrammes tombent au filtre de
continuite (>= 100 valeurs distinctes, atome <= 2 %). Une NLL est un flottant :
elle passe sans perte. La perplexite est donc le seul bloc a la fois dense,
continu, nombreux ET interpretable dans le cadre de la curation.

C'est aussi le design du papier arXiv 2509.23488, qui fait de la perplexite par
jeton ses features -- mais avec la moyenne conditionnelle pour cible. On garde
leurs features et on change de regime : la cible devient l'indice de queue.

Decoupage temporel, pour que X ne soit pas une fonction deterministe de la
passe qui produit Y :

    jetons     0..63      amorce, jamais scoree
             64..319      PREFIXE -> covariables X (256 jetons)
            320..831      SUFFIXE -> reponse Y     (512 jetons)

Y = exp(max_t NLL_t) sur le suffixe sous MODEL_Y.
X = profil de NLL sur le prefixe sous chacun de MODELS_X, plus les ecarts
    entre modeles -- l'analogue du gamma_Delta de la brique de curation.

usage: python code/py/extract_ppl_covariates.py DOCS.parquet OUTDIR [N] [GPU]
"""
from __future__ import annotations

import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")

DOCS = sys.argv[1] if len(sys.argv) > 1 else "results/wild/main_docs.parquet"
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else "results/wild/ppl"
NDOC = int(sys.argv[3]) if len(sys.argv) > 3 else 5000
GPU = sys.argv[4] if len(sys.argv) > 4 else "1"

os.environ["CUDA_VISIBLE_DEVICES"] = GPU
os.environ.setdefault("HF_HOME", "/people/pommeret/.cache/huggingface")
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")

import numpy as np
import pandas as pd
import torch
from transformers import AutoModelForCausalLM, AutoTokenizer

BURN, NPRE, NSUF = 64, 256, 512
NEED = BURN + NPRE + NSUF            # 832 jetons
BATCH = 8
MODEL_Y = "EleutherAI/pythia-410m"
MODELS_X = ["EleutherAI/pythia-70m", "EleutherAI/pythia-160m",
            "EleutherAI/pythia-410m"]
QLEV = np.round(np.linspace(0.02, 0.98, 49), 4)
NBLOCK = 16                          # blocs de 16 jetons sur le prefixe
NFBUCK = 20                          # tranches d'identifiant de jeton


def profile(nll: np.ndarray, tid: np.ndarray, tag: str) -> dict:
    """Toutes les statistiques d'un profil de NLL. Continues par construction."""
    f = {}
    q = np.quantile(nll, QLEV)
    for lev, v in zip(QLEV, q):
        f[f"{tag}_q{int(lev*100):02d}"] = float(v)
    f[f"{tag}_mean"] = float(nll.mean())
    f[f"{tag}_sd"] = float(nll.std())
    f[f"{tag}_max"] = float(nll.max())
    f[f"{tag}_min"] = float(nll.min())
    c = nll - nll.mean()
    s = nll.std() + 1e-12
    f[f"{tag}_skew"] = float((c ** 3).mean() / s ** 3)
    f[f"{tag}_kurt"] = float((c ** 4).mean() / s ** 4)
    # dynamique locale : moyennes et ecarts-types par bloc de 16 jetons
    B = nll.reshape(NBLOCK, -1)
    for b in range(NBLOCK):
        f[f"{tag}_bm{b:02d}"] = float(B[b].mean())
        f[f"{tag}_bs{b:02d}"] = float(B[b].std())
    # les identifiants BPE de GPT-NeoX sont grossierement ordonnes par
    # frequence : la tranche sert de proxy de rarete lexicale
    edges = np.linspace(0, 50304, NFBUCK + 1)
    idx = np.clip(np.digitize(tid, edges[1:-1]), 0, NFBUCK - 1)
    for b in range(NFBUCK):
        m = idx == b
        f[f"{tag}_fb{b:02d}"] = float(nll[m].mean()) if m.sum() >= 3 else np.nan
    return f


def main() -> None:
    os.makedirs(OUTDIR, exist_ok=True)
    df = pd.read_parquet(DOCS, columns=["doc_id", "text"])
    tok = AutoTokenizer.from_pretrained(MODEL_Y)

    # --- fenetrage commun a tous les modeles ---------------------------------
    t0 = time.time()
    ids, toks = [], []
    for s in range(0, len(df), 2000):
        part = df.iloc[s:s + 2000]
        enc = tok(part["text"].tolist(), truncation=True,
                  max_length=NEED)["input_ids"]
        for e, did in zip(enc, part["doc_id"].values):
            if len(e) >= NEED:
                toks.append(np.asarray(e[:NEED], dtype=np.int64))
                ids.append(int(did))
        if len(toks) >= NDOC:
            break
    toks, ids = toks[:NDOC], ids[:NDOC]
    n = len(toks)
    Tk = np.stack(toks)
    print(f"{n} fenetres de {NEED} jetons ({time.time()-t0:.0f} s)")

    def nll_all(model_name):
        """NLL par jeton pour tous les documents, positions 1..NEED-1."""
        m = AutoModelForCausalLM.from_pretrained(
            model_name, torch_dtype=torch.float32).cuda().eval()
        out = np.empty((n, NEED - 1), dtype=np.float32)
        t = time.time()
        for s in range(0, n, BATCH):
            x = torch.tensor(Tk[s:s + BATCH], dtype=torch.long, device="cuda")
            with torch.no_grad():
                lg = m(x).logits.float()
                lp = torch.log_softmax(lg[:, :-1], dim=-1)
                v = -lp.gather(2, x[:, 1:, None]).squeeze(2)
            out[s:s + BATCH] = v.cpu().numpy()
        del m; torch.cuda.empty_cache()
        print(f"  {model_name} : {time.time()-t:.0f} s")
        return out

    # --- reponse : max sur le suffixe sous MODEL_Y ---------------------------
    NY = nll_all(MODEL_Y)
    suf = NY[:, BURN + NPRE - 1:]
    assert suf.shape[1] == NSUF, suf.shape
    Ylog = suf.max(axis=1)

    # --- covariables : profils de prefixe -----------------------------------
    feats, pre_store = [{} for _ in range(n)], {}
    tid_pre = Tk[:, BURN:BURN + NPRE]   # aligne sur pre = NX[:, BURN-1:...]
    for mn in MODELS_X:
        NX = NY if mn == MODEL_Y else nll_all(mn)
        pre = NX[:, BURN - 1:BURN + NPRE - 1]
        assert pre.shape[1] == NPRE, pre.shape
        pre_store[mn] = pre
        tag = "p" + mn.split("-")[-1]
        for i in range(n):
            feats[i].update(profile(pre[i].astype(np.float64), tid_pre[i], tag))

    # --- ecarts entre modeles : l'analogue de gamma_Delta --------------------
    for a, b in [("EleutherAI/pythia-70m", "EleutherAI/pythia-410m"),
                 ("EleutherAI/pythia-160m", "EleutherAI/pythia-410m")]:
        d = pre_store[a] - pre_store[b]
        tag = f"d{a.split('-')[-1]}_{b.split('-')[-1]}"
        for i in range(n):
            feats[i].update(profile(d[i].astype(np.float64), tid_pre[i], tag))

    F = pd.DataFrame(feats)
    F.insert(0, "doc_id", ids)
    F["Ylog"] = Ylog
    F["suffix_mean"] = suf.mean(axis=1)
    bad = F.isna().any(axis=0)
    if bad.any():
        print(f"  colonnes a NaN ecartees : {int(bad.sum())}")
        F = F.loc[:, ~bad]
    F.to_parquet(os.path.join(OUTDIR, "ppl_cov.parquet"), index=False)
    ncov = F.shape[1] - 3
    print(f"ECRIT {OUTDIR}/ppl_cov.parquet : {n} x {ncov} covariables")
    print(f"  Ylog : med {np.median(Ylog):.3f}, q99 {np.quantile(Ylog,.99):.3f}, "
          f"max {Ylog.max():.3f}")


if __name__ == "__main__":
    main()
