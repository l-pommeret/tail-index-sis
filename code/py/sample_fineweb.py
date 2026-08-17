"""Tirage de documents FineWeb, avec plancher de longueur en jetons.

Le plancher est indispensable : extract_nll.py score EXACTEMENT T = 512 jetons
apres 64 d'amorce, donc un document de moins de 576 jetons ne peut pas entrer
sans reintroduire la longueur comme nuisance d'echelle. token_count du parquet
vient d'un autre tokenizer que Pythia, on prend donc une marge.

usage: python code/py/sample_fineweb.py OUT.parquet [N] [SEED]
"""
import os, sys, warnings; warnings.filterwarnings("ignore")
os.environ.setdefault("HF_HOME", "/people/pommeret/.cache/huggingface")
import numpy as np, pandas as pd, pyarrow.parquet as pq
from huggingface_hub import hf_hub_download

OUT = sys.argv[1] if len(sys.argv) > 1 else "results/wild/docs.parquet"
N = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 20260817
MIN_TOK = 700          # marge sur les 576 jetons Pythia requis

p = hf_hub_download("HuggingFaceFW/fineweb", "sample/100BT/000_00000.parquet",
                    repo_type="dataset")
f = pq.ParquetFile(p)
rng = np.random.default_rng(SEED)
keep = []
for b in f.iter_batches(batch_size=50000,
                        columns=["text", "id", "url", "dump", "token_count",
                                 "language_score"]):
    d = b.to_pandas()
    d = d[d.token_count >= MIN_TOK]
    keep.append(d)
    if sum(len(x) for x in keep) >= 6 * N:
        break
d = pd.concat(keep, ignore_index=True)
idx = rng.choice(len(d), size=min(N, len(d)), replace=False)
d = d.iloc[np.sort(idx)].reset_index(drop=True)
d["doc_id"] = np.arange(len(d))
os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
d.to_parquet(OUT, index=False)
print(f"ECRIT {OUT} : {len(d)} documents (plancher {MIN_TOK} jetons)")
print(f"  token_count : med {d.token_count.median():.0f}, "
      f"q10 {d.token_count.quantile(.1):.0f}, q90 {d.token_count.quantile(.9):.0f}")
