"""Etats internes du modele comme covariables, reponse sur le futur.

Question : quelles directions du flux residuel annoncent que le modele est sur
le point d'echouer catastrophiquement ?

Le point 5 du protocole interdit de tirer X de la passe forward qui produit Y :
si X determine Y, la loi conditionnelle degenere. La parade est un decoupage
passe / futur dans le meme document,

    X_i = activations moyennees sur les tokens du PREFIXE  (positions B..B+T0)
    Y_i = exp( max_t S_t ) sur le SUFFIXE                  (positions B+T0..B+T0+T)

X est mesurable par rapport au passe, Y par rapport au futur : c'est un
probleme de prediction, pas une tautologie.

Huit couches sont conservees (3, 6, 9, 12, 15, 18, 21, 24 sur 24) plutot que
trois : le flux residuel evolue lentement, donc des couches voisines sont quasi
colineaires et n'ajoutent pas de dimensions effectives -- le meme phenomene que
sur les grilles de reglage, ou densifier n'augmentait pas le nombre effectif de
configurations independantes. En garder huit bien separees permet de tester
plusieurs configurations sans reextraire.

usage: python code/py/extract_internal.py OUT.parquet [N_DOCS]
"""
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import torch
from datasets import load_dataset
from transformers import AutoModelForCausalLM, AutoTokenizer

OUT = sys.argv[1] if len(sys.argv) > 1 else "results/llm/internal.parquet"
N_DOCS = int(sys.argv[2]) if len(sys.argv) > 2 else 20000

MODEL = "EleutherAI/pythia-410m-deduped"
DUMP = "CC-MAIN-2025-18"
BURN_IN = 64          # tokens de tete ecartes
T0 = 256              # longueur du prefixe, qui porte X
T = 512               # longueur du suffixe, qui porte Y
CTX = 2048
LAYERS = [3, 6, 9, 12, 15, 18, 21, 24]


@torch.no_grad()
def one_doc(model, ids, dev, special):
    """Retourne (activations du prefixe par couche, max de surprise du suffixe)."""
    out = model(ids.to(dev), output_hidden_states=True)
    lg = out.logits.float()
    lp = torch.log_softmax(lg[0, :-1], dim=-1)
    S = -lp.gather(1, ids[0, 1:, None].to(dev)).squeeze(1)
    tgt = ids[0, 1:]
    ok = torch.tensor([int(t) not in special for t in tgt], device=dev)

    pre = slice(BURN_IN, BURN_IN + T0)
    suf = slice(BURN_IN + T0, BURN_IN + T0 + T)
    s_suf = S[suf][ok[suf]]
    if s_suf.numel() < T // 2:
        return None, None
    # hidden_states[0] est la sortie de l'embedding, [i] celle de la couche i
    acts = [out.hidden_states[l][0, pre].mean(0).float().cpu().numpy()
            for l in LAYERS]
    return np.concatenate(acts), float(s_suf.max())


def _save(ids_l, acts_l, y_l, path):
    A = np.vstack(acts_l).astype(np.float32)
    per = A.shape[1] // len(LAYERS)
    cols = [f"h{l:02d}_{d:04d}" for l in LAYERS for d in range(per)]
    df = pd.DataFrame(A, columns=cols)
    df.insert(0, "doc_id", ids_l)
    df.insert(1, "max_S_suffix", y_l)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_parquet(path, index=False)
    return df


def main():
    tok = AutoTokenizer.from_pretrained(MODEL)
    special = set(tok.all_special_ids)
    model = AutoModelForCausalLM.from_pretrained(MODEL, dtype=torch.float32).cuda(0).eval()
    ds = load_dataset("HuggingFaceFW/fineweb", DUMP, split="train", streaming=True)

    need = BURN_IN + T0 + T + 1
    ids_l, acts_l, y_l, t0, seen = [], [], [], time.time(), 0
    for ex in ds:
        if len(y_l) >= N_DOCS:
            break
        seen += 1
        ids = tok(ex["text"], return_tensors="pt").input_ids
        if ids.shape[1] < need:
            continue
        a, y = one_doc(model, ids[:, :CTX], 0, special)
        if a is None:
            continue
        ids_l.append(ex.get("id", "")); acts_l.append(a); y_l.append(y)
        if len(y_l) % 1000 == 0:
            el = time.time() - t0
            print(f"  {len(y_l)}/{N_DOCS} ({seen} lus) {el:.0f}s "
                  f"[{len(y_l)/el:.1f} doc/s]", flush=True)
        if len(y_l) % 10000 == 0:      # point de reprise
            _save(ids_l, acts_l, y_l, OUT + ".part")

    df = _save(ids_l, acts_l, y_l, OUT)
    nact = df.shape[1] - 2
    print(f"\nECRIT {OUT} : {len(df)} documents, {nact} activations "
          f"({len(LAYERS)} couches x {nact // len(LAYERS)})")
    print(f"documents lus : {seen}")
    print(f"max_S sur le suffixe : med {np.median(y_l):.2f}, "
          f"q99 {np.quantile(y_l, .99):.2f}, max {max(y_l):.2f} nats")


if __name__ == "__main__":
    main()
