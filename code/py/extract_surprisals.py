"""Extraction des surprises par token pour l'application LLM.

Protocole (points 3 et 5) :
  - modeles geles, fp32 de bout en bout, logits castes en fp32 avant softmax ;
  - batch = 1, donc la numerique des kernels ne depend pas du remplissage ;
  - burn-in de B tokens ecartes en tete, tokens speciaux exclus ;
  - deux modeles de la meme famille, donc MEME tokenizer : les surprises sont
    alignees token par token et la difference S_petit - S_grand a un sens ;
  - covariables calculees a partir du texte brut uniquement, jamais de la passe
    forward qui produit la reponse (response-blind).

Deux reponses sont enregistrees par document :
  Y      = exp(max_t S_t)            indice de queue du modele evalue
  Y_delta= exp(max_t (S_petit - S_grand))   rapport de vraisemblance maximal

Deux decoupages :
  fixe  : exactement T tokens notes par document (analyse primaire, aucun
          confondant de longueur)
  libre : le document entier jusqu'a la limite de contexte (analyse secondaire,
          la longueur devient la nuisance d'echelle)

usage: python code/py/extract_surprisals.py OUT.parquet [N_DOCS] [MODE]
       MODE = fixe | libre
"""
import gzip
import json
import math
import os
import re
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import torch
from datasets import load_dataset
from transformers import AutoModelForCausalLM, AutoTokenizer

OUT = sys.argv[1] if len(sys.argv) > 1 else "results/llm/pilot.parquet"
N_DOCS = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
MODE = sys.argv[3] if len(sys.argv) > 3 else "fixe"

SMALL = "EleutherAI/pythia-410m-deduped"
LARGE = "EleutherAI/pythia-1.4b-deduped"
DUMP = "CC-MAIN-2025-18"          # crawl tres posterieur au Pile (collecte 2020)
BURN_IN = 64                       # tokens de tete ecartes
T_FIXED = 512                      # tokens notes en mode fixe
CTX = 2048                         # limite de contexte des Pythia
TOPK_KEEP = 8                      # plus grandes surprises conservees par doc

# ------------------------------------------------------------------ covariables
WORD = re.compile(r"[A-Za-z]+")
SENT = re.compile(r"[.!?]+(?:\s|$)")


def surface_covariates(txt: str) -> dict:
    """Statistiques de surface, calculees sans jamais toucher au modele.

    Le bloc 'longueur' est deliberement redondant : caracteres, octets, mots,
    phrases, paragraphes, espaces, taille compressee. La theorie des valeurs
    extremes donne P(max > y) ~ theta T c y^(-1/gamma), donc la longueur entre
    dans l'echelle et jamais dans l'indice : ces coordonnees forment un amas de
    proxys d'une nuisance d'echelle dont le statut est etabli analytiquement.
    """
    b = txt.encode("utf-8", "ignore")
    words = WORD.findall(txt)
    lw = [w.lower() for w in words]
    n_ch = max(len(txt), 1)
    n_w = max(len(words), 1)
    uniq = set(lw)
    counts = {}
    for w in lw:
        counts[w] = counts.get(w, 0) + 1
    hapax = sum(1 for v in counts.values() if v == 1)
    comp = len(gzip.compress(b, 6))
    d = {
        # --- amas longueur (nuisance d'echelle, verite analytique) ---
        "len_chars": len(txt),
        "len_bytes": len(b),
        "len_words": len(words),
        "len_sents": len(SENT.findall(txt)) + 1,
        "len_paras": txt.count("\n\n") + 1,
        "len_spaces": txt.count(" "),
        "len_lines": txt.count("\n") + 1,
        "len_gzip": comp,
        # --- typographie et composition ---
        "rate_digit": sum(c.isdigit() for c in txt) / n_ch,
        "rate_upper": sum(c.isupper() for c in txt) / n_ch,
        "rate_punct": sum(c in ".,;:!?'\"()[]{}-—/\\|" for c in txt) / n_ch,
        "rate_nonascii": sum(ord(c) > 127 for c in txt) / n_ch,
        "rate_space": txt.count(" ") / n_ch,
        "rate_newline": txt.count("\n") / n_ch,
        "rate_alpha": sum(c.isalpha() for c in txt) / n_ch,
        # --- lexique ---
        "ttr": len(uniq) / n_w,
        "hapax_rate": hapax / n_w,
        "mean_word_len": float(np.mean([len(w) for w in words])) if words else 0.0,
        "max_word_len": max((len(w) for w in words), default=0),
        "mean_sent_len": n_w / max(len(SENT.findall(txt)) + 1, 1),
        # --- compression et entropie de surface ---
        "gzip_ratio": comp / max(len(b), 1),
        "byte_entropy": float(
            -sum(
                (c / len(b)) * math.log2(c / len(b))
                for c in np.bincount(np.frombuffer(b, dtype=np.uint8), minlength=256)
                if c
            )
        )
        if b
        else 0.0,
    }
    return d


# ------------------------------------------------------------------- modeles --
def load_models():
    tok = AutoTokenizer.from_pretrained(SMALL)
    tok_l = AutoTokenizer.from_pretrained(LARGE)
    assert tok.get_vocab() == tok_l.get_vocab(), "tokenizers differents"
    ms = AutoModelForCausalLM.from_pretrained(SMALL, dtype=torch.float32).cuda(0).eval()
    ml = AutoModelForCausalLM.from_pretrained(LARGE, dtype=torch.float32).cuda(1).eval()
    return tok, ms, ml


@torch.no_grad()
def surprisals(model, ids, dev):
    """Surprises en nats pour les positions 1..T-1, en fp32."""
    lg = model(ids.to(dev)).logits.float()
    lp = torch.log_softmax(lg[0, :-1], dim=-1)
    return -lp.gather(1, ids[0, 1:, None].to(dev)).squeeze(1).float().cpu()


def main():
    tok, ms, ml = load_models()
    special = set(tok.all_special_ids)
    ds = load_dataset("HuggingFaceFW/fineweb", DUMP, split="train", streaming=True)

    rows, t0, seen = [], time.time(), 0
    need = BURN_IN + (T_FIXED if MODE == "fixe" else 128) + 1
    for ex in ds:
        if len(rows) >= N_DOCS:
            break
        seen += 1
        txt = ex["text"]
        ids = tok(txt, return_tensors="pt").input_ids
        if ids.shape[1] < need:
            continue
        ids = ids[:, :CTX]
        Ss = surprisals(ms, ids, 0)
        Sl = surprisals(ml, ids, 1)
        # positions notees : apres burn-in, hors tokens speciaux
        tgt = ids[0, 1:]
        keep = torch.tensor([i >= BURN_IN and int(t) not in special
                             for i, t in enumerate(tgt)])
        if MODE == "fixe":
            idx = torch.nonzero(keep).squeeze(1)[:T_FIXED]
            if idx.numel() < T_FIXED:
                continue
        else:
            idx = torch.nonzero(keep).squeeze(1)
            if idx.numel() < 64:
                continue
        s_small, s_large = Ss[idx], Sl[idx]
        delta = s_small - s_large
        top = torch.topk(s_small, min(TOPK_KEEP, s_small.numel())).values
        r = {
            "doc_id": ex.get("id", ""),
            "url": ex.get("url", ""),
            "dump": ex.get("dump", ""),
            "date": str(ex.get("date", "")),
            "n_tok_doc": int(ids.shape[1]),
            "T_scored": int(idx.numel()),
            "max_S_small": float(s_small.max()),
            "max_S_large": float(s_large.max()),
            "max_delta": float(delta.max()),
            "mean_S_small": float(s_small.mean()),
            "mean_S_large": float(s_large.mean()),
            "argmax_tok": tok.decode([int(tgt[idx[int(s_small.argmax())]])]),
            "top_S": json.dumps([round(float(v), 4) for v in top]),
        }
        r.update(surface_covariates(txt))
        rows.append(r)
        if len(rows) % 200 == 0:
            el = time.time() - t0
            print(f"{len(rows)}/{N_DOCS} docs ({seen} lus) {el:.0f}s "
                  f"[{len(rows)/el:.1f} doc/s]", flush=True)

    import pandas as pd
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    df.to_parquet(OUT, index=False)
    print(f"\nECRIT {OUT} : {len(df)} documents, {df.shape[1]} colonnes")
    print(f"documents lus pour en retenir {len(df)} : {seen}")
    print(f"max_S_small : med {df.max_S_small.median():.2f}, "
          f"q99 {df.max_S_small.quantile(.99):.2f}, max {df.max_S_small.max():.2f} nats")
    print(f"max_delta   : med {df.max_delta.median():.2f}, "
          f"max {df.max_delta.max():.2f} nats")


if __name__ == "__main__":
    main()
