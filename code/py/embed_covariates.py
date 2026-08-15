"""Covariables d'embedding, depuis un encodeur separe et gele.

Point 5 du protocole : les covariables ne doivent pas provenir de la passe
forward qui produit la reponse. Si X etait une fonction deterministe des memes
logits que Y, la loi conditionnelle pourrait degenerer. On utilise donc un
encodeur d'une autre famille, gele, et on ne lui demande que des embeddings de
surface -- jamais de log-probabilites.

Le flux et le filtre de selection sont identiques a extract_surprisals.py, et la
lecture du corpus est deterministe, donc les memes documents reviennent dans le
meme ordre ; l'identifiant est conserve pour verifier l'appariement.

usage: python code/py/embed_covariates.py OUT.parquet [N_DOCS] [N_PCA]
"""
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import torch
from datasets import load_dataset
from transformers import AutoModel, AutoTokenizer

OUT = sys.argv[1] if len(sys.argv) > 1 else "results/llm/main_embed.parquet"
N_DOCS = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
N_PCA = int(sys.argv[3]) if len(sys.argv) > 3 else 128

ENCODER = "roberta-large"          # famille distincte des Pythia qui produisent Y
GEN_TOK = "EleutherAI/pythia-410m-deduped"   # pour reproduire le filtre de selection
DUMP = "CC-MAIN-2025-18"
BURN_IN, T_FIXED, CTX = 64, 512, 2048
BATCH = 16


@torch.no_grad()
def embed(model, tok, texts, dev):
    b = tok(texts, return_tensors="pt", padding=True, truncation=True,
            max_length=512).to(dev)
    h = model(**b).last_hidden_state
    mask = b["attention_mask"].unsqueeze(-1).float()
    return ((h * mask).sum(1) / mask.sum(1).clamp(min=1)).float().cpu().numpy()


def main():
    gen_tok = AutoTokenizer.from_pretrained(GEN_TOK)
    tok = AutoTokenizer.from_pretrained(ENCODER)
    model = AutoModel.from_pretrained(ENCODER, dtype=torch.float32).cuda(0).eval()
    need = BURN_IN + T_FIXED + 1

    ds = load_dataset("HuggingFaceFW/fineweb", DUMP, split="train", streaming=True)
    ids, texts, embs = [], [], []
    buf_id, buf_tx = [], []
    for ex in ds:
        if len(ids) >= N_DOCS:
            break
        n_tok = gen_tok(ex["text"], return_tensors="pt").input_ids.shape[1]
        if n_tok < need:
            continue
        buf_id.append(ex.get("id", "")); buf_tx.append(ex["text"])
        if len(buf_tx) == BATCH:
            embs.append(embed(model, tok, buf_tx, 0))
            ids.extend(buf_id); buf_id, buf_tx = [], []
            if len(ids) % 2000 == 0:
                print(f"  {len(ids)}/{N_DOCS}", flush=True)
    if buf_tx:
        embs.append(embed(model, tok, buf_tx, 0)); ids.extend(buf_id)

    E = np.vstack(embs)[:len(ids)]
    print(f"embeddings bruts : {E.shape}")

    # ACP : on garde les composantes qui portent le signal, pas les 1024 brutes
    Ec = E - E.mean(0, keepdims=True)
    U, S, Vt = np.linalg.svd(Ec, full_matrices=False)
    k = min(N_PCA, Vt.shape[0])
    P = Ec @ Vt[:k].T
    var = (S ** 2) / (S ** 2).sum()
    print(f"ACP {k} composantes, variance expliquee {var[:k].sum():.3f}")

    df = pd.DataFrame(P, columns=[f"emb_{i:03d}" for i in range(k)])
    df.insert(0, "doc_id", ids)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    df.to_parquet(OUT, index=False)
    print(f"ECRIT {OUT} : {len(df)} documents, {k} composantes")


if __name__ == "__main__":
    main()
