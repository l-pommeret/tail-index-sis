#!/bin/sh
# Enchaine les etapes de l'application LLM apres l'extraction principale :
#   1. embeddings d'un encodeur separe et gele (roberta-large), ACP
#   2. jointure des covariables de surface et d'embedding, export CSV
#   3. criblage, calibration par permutation, frequences de selection
set -u
cd /people/pommeret/tail-index-sis
PY=/people/pommeret/miniconda3/bin/python
export PATH=/people/pommeret/tail-index-sis/.renv/bin:$PATH

until [ -f results/llm/main_fixe.parquet ]; do sleep 60; done
echo "$(date +%H:%M:%S) extraction principale terminee"

$PY code/py/embed_covariates.py results/llm/main_embed.parquet 20000 256 \
    > results/llm/main_embed.log 2>&1
echo "$(date +%H:%M:%S) embeddings termines"

$PY - <<'PYEOF' > results/llm/merge.log 2>&1
import pandas as pd
a = pd.read_parquet("results/llm/main_fixe.parquet")
b = pd.read_parquet("results/llm/main_embed.parquet")
m = a.merge(b, on="doc_id", how="inner")
print(f"surface {a.shape}, embed {b.shape}, jointure {m.shape}")
print(f"documents perdus a la jointure : {len(a) - len(m)}")
m.to_csv("results/llm/main_fixe.csv", index=False)
PYEOF
echo "$(date +%H:%M:%S) jointure faite"

Rscript code/R/screen_llm.R results/llm/main_fixe.csv results/llm/screen 92 1000 \
    > results/llm/screen.log 2>&1
echo "$(date +%H:%M:%S) criblage termine"
