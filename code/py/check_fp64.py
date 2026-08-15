"""Verification fp64 des plus grandes surprises (point 3 du protocole).

Les plus grandes surprises sont les observations que la precision reduite
degrade le plus, et ce sont precisement elles qui pilotent Hill. On reprend les
K documents de plus grand max_t S, on recalcule leurs surprises avec le modele
et le softmax en float64, et on verifie que le classement est inchange.

Les documents sont retrouves en rejouant le flux du corpus et en filtrant sur
l'identifiant : la lecture est deterministe, donc les memes documents
reviennent dans le meme ordre.

usage: python code/py/check_fp64.py PILOT.parquet [K]
"""
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import torch
from datasets import load_dataset
from scipy.stats import kendalltau, spearmanr
from transformers import AutoModelForCausalLM, AutoTokenizer

PATH = sys.argv[1] if len(sys.argv) > 1 else "results/llm/pilot_fixe.parquet"
K = int(sys.argv[2]) if len(sys.argv) > 2 else 200

SMALL = "EleutherAI/pythia-410m-deduped"
DUMP = "CC-MAIN-2025-18"
BURN_IN, T_FIXED, CTX = 64, 512, 2048


@torch.no_grad()
def max_surprisal(model, ids, dev, special, dtype):
    lg = model(ids.to(dev)).logits.to(dtype)
    lp = torch.log_softmax(lg[0, :-1], dim=-1)
    S = -lp.gather(1, ids[0, 1:, None].to(dev)).squeeze(1)
    tgt = ids[0, 1:]
    keep = torch.tensor([i >= BURN_IN and int(t) not in special
                         for i, t in enumerate(tgt)])
    idx = torch.nonzero(keep).squeeze(1)[:T_FIXED]
    return float(S[idx.to(dev)].max())


def main():
    df = pd.read_parquet(PATH)
    top = df.nlargest(K, "max_S_small")
    want = dict(zip(top.doc_id, top.max_S_small))
    print(f"{K} documents de plus grand max_t S, de {top.max_S_small.min():.2f} "
          f"a {top.max_S_small.max():.2f} nats")

    tok = AutoTokenizer.from_pretrained(SMALL)
    special = set(tok.all_special_ids)
    m32 = AutoModelForCausalLM.from_pretrained(SMALL, dtype=torch.float32).cuda(0).eval()
    m64 = AutoModelForCausalLM.from_pretrained(SMALL, dtype=torch.float64).cuda(1).eval()
    print("modeles charges en fp32 et fp64")

    ds = load_dataset("HuggingFaceFW/fineweb", DUMP, split="train", streaming=True)
    rows, seen = [], 0
    for ex in ds:
        seen += 1
        if ex.get("id") not in want:
            if len(rows) >= len(want):
                break
            continue
        ids = tok(ex["text"], return_tensors="pt").input_ids[:, :CTX]
        if ids.shape[1] < BURN_IN + T_FIXED + 1:
            continue
        s32 = max_surprisal(m32, ids, 0, special, torch.float32)
        s64 = max_surprisal(m64, ids, 1, special, torch.float64)
        rows.append(dict(doc_id=ex["id"], stored=want[ex["id"]], fp32=s32, fp64=s64))
        if len(rows) % 50 == 0:
            print(f"  {len(rows)}/{len(want)} recalcules ({seen} lus)", flush=True)
        if len(rows) >= len(want):
            break

    r = pd.DataFrame(rows)
    r["ecart"] = r.fp64 - r.fp32
    tau = kendalltau(r.fp32, r.fp64).statistic
    rho = spearmanr(r.fp32, r.fp64).statistic
    inv32 = np.argsort(-r.fp32.values)
    inv64 = np.argsort(-r.fp64.values)
    print(f"\n{len(r)} documents recalcules")
    print(f"reproduction du stockage fp32 : ecart max {np.abs(r.stored - r.fp32).max():.2e}")
    print(f"ecart fp64 - fp32 : med {r.ecart.median():+.3e}, "
          f"max abs {r.ecart.abs().max():.3e} nats")
    print(f"correlation de rang fp32 vs fp64 : Kendall {tau:.6f}, Spearman {rho:.6f}")
    print(f"rangs identiques : {int((inv32 == inv64).sum())}/{len(r)}")
    print(f"meme document en tete : {inv32[0] == inv64[0]}")
    top20 = set(inv32[:20]) & set(inv64[:20])
    print(f"intersection des 20 premiers : {len(top20)}/20")
    r.to_csv(PATH.replace(".parquet", "_fp64.csv"), index=False)


if __name__ == "__main__":
    main()
