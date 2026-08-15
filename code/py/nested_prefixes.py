"""Effet pur de la longueur : prefixes emboites d'un meme document.

La comparaison entre strates de longueur melange deux choses : l'effet de T et
le fait que les documents longs different des courts par leur contenu. En
scorant des PREFIXES EMBOITES du meme document, la composition est tenue
exactement fixe et la pente mesure le seul effet de T.

La theorie des valeurs extremes donne, pour un document x,

    E[M_T | x] = gamma(x) log T + O(1),

donc la pente de M_T contre log T, estimee DANS chaque document puis moyennee,
estime gamma sans aucun confondant de composition. C'est la verification la plus
propre du statut de la longueur comme nuisance d'echelle.

usage: python code/py/nested_prefixes.py OUT.csv [N_DOCS]
"""
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import torch
from datasets import load_dataset
from transformers import AutoModelForCausalLM, AutoTokenizer

OUT = sys.argv[1] if len(sys.argv) > 1 else "results/llm/nested.csv"
N_DOCS = int(sys.argv[2]) if len(sys.argv) > 2 else 600

SMALL = "EleutherAI/pythia-410m-deduped"
DUMP = "CC-MAIN-2025-18"
BURN_IN = 64
TS = [128, 256, 512, 1024]        # longueurs emboitees
CTX = 2048


@torch.no_grad()
def prefix_maxima(model, ids, dev, special):
    """max_t S_t sur les TS premiers tokens notes, en une seule passe."""
    lg = model(ids.to(dev)).logits.float()
    lp = torch.log_softmax(lg[0, :-1], dim=-1)
    S = -lp.gather(1, ids[0, 1:, None].to(dev)).squeeze(1).cpu()
    tgt = ids[0, 1:]
    keep = torch.tensor([i >= BURN_IN and int(t) not in special
                         for i, t in enumerate(tgt)])
    idx = torch.nonzero(keep).squeeze(1)
    return [float(S[idx[:T]].max()) for T in TS]


def main():
    tok = AutoTokenizer.from_pretrained(SMALL)
    special = set(tok.all_special_ids)
    m = AutoModelForCausalLM.from_pretrained(SMALL, dtype=torch.float32).cuda(0).eval()
    ds = load_dataset("HuggingFaceFW/fineweb", DUMP, split="train", streaming=True)

    need = BURN_IN + max(TS) + 1
    rows = []
    for ex in ds:
        if len(rows) >= N_DOCS:
            break
        ids = tok(ex["text"], return_tensors="pt").input_ids
        if ids.shape[1] < need:
            continue
        ids = ids[:, :CTX]
        ms = prefix_maxima(m, ids, 0, special)
        rows.append(dict(doc_id=ex.get("id", ""),
                         **{f"M_{T}": v for T, v in zip(TS, ms)}))
        if len(rows) % 100 == 0:
            print(f"  {len(rows)}/{N_DOCS}", flush=True)

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    df.to_csv(OUT, index=False)

    lt = np.log(np.array(TS, dtype=float))
    lt_c = lt - lt.mean()
    denom = float((lt_c ** 2).sum())
    Mmat = df[[f"M_{T}" for T in TS]].values
    slopes = ((Mmat - Mmat.mean(axis=1, keepdims=True)) @ lt_c) / denom

    print(f"\n{len(df)} documents, prefixes T = {TS}")
    print("\nmoyenne de M par longueur (composition tenue fixe) :")
    for T, v in zip(TS, Mmat.mean(axis=0)):
        print(f"  T={T:5d}  log T={np.log(T):.3f}  E[M] = {v:.3f} nats")
    dmean = Mmat.mean(axis=0)
    A = np.vstack([lt, np.ones_like(lt)]).T
    sl, ic = np.linalg.lstsq(A, dmean, rcond=None)[0]
    resid = dmean - A @ [sl, ic]
    print(f"\npente sur les moyennes : {sl:.4f} nats/nat "
          f"(residu max {np.abs(resid).max():.4f})")
    print(f"pente intra-document : moyenne {slopes.mean():.4f}, "
          f"mediane {np.median(slopes):.4f}, "
          f"erreur type {slopes.std(ddof=1)/np.sqrt(len(slopes)):.4f}")
    lo = slopes.mean() - 1.96 * slopes.std(ddof=1) / np.sqrt(len(slopes))
    hi = slopes.mean() + 1.96 * slopes.std(ddof=1) / np.sqrt(len(slopes))
    print(f"IC95 de la pente intra-document : [{lo:.4f}, {hi:.4f}]")
    print("\nla theorie predit pente = gamma ; comparer a Hill ~ 1.30-1.38")
    print("croissance monotone par document : "
          f"{float(np.mean(np.all(np.diff(Mmat, axis=1) >= 0, axis=1))):.3f} "
          "(M_T est un maximum sur un ensemble croissant, donc doit etre non decroissant)")


if __name__ == "__main__":
    main()
