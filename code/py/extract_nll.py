"""Surprisal par jeton sous Pythia, et reponse Y = exp(max NLL) par document.

Pourquoi exp(max) et pas la perplexite moyenne. Si la surprisal par jeton
S_t = -log p_theta(x_t | x_<t) a une queue superieure de type exponentiel,
P(S > x) ~ C e^{-lambda x}, alors Y = e^S verifie P(Y > y) ~ C y^{-lambda} :
Pareto EXACTEMENT, d'indice gamma = 1/lambda. Le maximum sur T jetons conserve
l'indice (le max d'un echantillon a variation reguliere l'est aussi) et donne
E[max_t S_t] = gamma log T + O(1). L'indice de queue de Y se lit donc comme le
taux de decroissance de la distribution de surprise du modele, et gamma(x)
grand signifie "ce profil de covariables admet des jetons catastrophiquement
surprenants" -- le critere de curation.

Deux precautions de mesure :

  amorce   les premiers jetons d'un document ont une surprisal artificiellement
           haute (le modele n'a pas de contexte). On saute BURN jetons.
  T fixe   la longueur est une nuisance d'echelle : E[max] croit en gamma log T.
           On score EXACTEMENT T jetons par document, ce qui la neutralise par
           construction. Les documents trop courts sont ecartes en amont.

usage: python code/py/extract_nll.py IN.parquet OUT.parquet [N] [MODEL] [GPU]
"""
from __future__ import annotations

import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")

IN = sys.argv[1] if len(sys.argv) > 1 else "results/wild/docs.parquet"
OUT = sys.argv[2] if len(sys.argv) > 2 else "results/wild/nll.parquet"
NDOC = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
MODEL = sys.argv[4] if len(sys.argv) > 4 else "EleutherAI/pythia-410m"
GPU = sys.argv[5] if len(sys.argv) > 5 else "1"

os.environ["CUDA_VISIBLE_DEVICES"] = GPU
os.environ.setdefault("HF_HOME", "/people/pommeret/.cache/huggingface")

import numpy as np
import pandas as pd
import torch
from transformers import AutoModelForCausalLM, AutoTokenizer

BURN = 64          # jetons d'amorce non scores
T = 512            # jetons scores, identique pour tous les documents
NEED = BURN + T    # longueur minimale en jetons
BATCH = 16


def main() -> None:
    df = pd.read_parquet(IN)
    print(f"{len(df)} documents en entree, cible {NDOC}")

    tok = AutoTokenizer.from_pretrained(MODEL)
    model = AutoModelForCausalLM.from_pretrained(
        MODEL, torch_dtype=torch.float32).cuda().eval()
    print(f"{MODEL} charge sur GPU {GPU}, fp32")

    ids_keep, rows, buf, meta = [], [], [], []
    t0 = time.time()

    def flush():
        if not buf:
            return
        x = torch.tensor(np.stack(buf), dtype=torch.long, device="cuda")
        with torch.no_grad():
            logits = model(x).logits.float()
            # NLL du jeton t predit depuis <t : on aligne logits[:, :-1] sur x[:, 1:]
            lp = torch.log_softmax(logits[:, :-1], dim=-1)
            nll = -lp.gather(2, x[:, 1:, None]).squeeze(2)   # (B, NEED-1)
            s = nll[:, BURN - 1:]                            # exactement T valeurs
        s = s.cpu().numpy()
        assert s.shape[1] == T, (s.shape, T)
        for k, m in enumerate(meta):
            v = s[k]
            rows.append({**m,
                         "max_nll": float(v.max()),
                         "mean_nll": float(v.mean()),
                         "sd_nll": float(v.std()),
                         "q99_nll": float(np.quantile(v, 0.99)),
                         "argmax_pos": int(v.argmax()),
                         "top2_nll": float(np.sort(v)[-2])})
        buf.clear(); meta.clear()

    for i, r in enumerate(df.itertuples()):
        e = tok(r.text, truncation=True, max_length=NEED)["input_ids"]
        if len(e) < NEED:
            continue
        buf.append(np.asarray(e[:NEED], dtype=np.int64))
        meta.append({"doc_id": int(getattr(r, "doc_id", i))})
        ids_keep.append(int(getattr(r, "doc_id", i)))
        if len(buf) == BATCH:
            flush()
        if len(rows) >= NDOC:
            break
        if len(rows) and len(rows) % 500 == 0 and not buf:
            el = time.time() - t0
            print(f"  {len(rows)} docs, {el:.0f} s, {len(rows)/el:.1f} doc/s",
                  flush=True)
    flush()

    out = pd.DataFrame(rows)
    out["Y"] = np.exp(out["max_nll"])
    os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
    out.to_parquet(OUT, index=False)
    el = time.time() - t0
    print(f"ECRIT {OUT} : {len(out)} documents en {el:.0f} s "
          f"({len(out)/el:.1f} doc/s)")
    print(f"  max_nll : med {out.max_nll.median():.3f}, "
          f"q99 {out.max_nll.quantile(.99):.3f}, max {out.max_nll.max():.3f}")
    print(f"  T = {T} jetons scores, amorce {BURN}, tous les documents identiques")


if __name__ == "__main__":
    main()
