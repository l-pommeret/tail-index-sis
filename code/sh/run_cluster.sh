#!/bin/sh
# Criblage des cinq bras a nuisance d'echelle injectee.
set -u
cd /people/pommeret/tail-index-sis
export PATH=/people/pommeret/tail-index-sis/.renv/bin:$PATH
PY=/people/pommeret/miniconda3/bin/python
for K in k000 k005 k010 k020 k040; do
  $PY -c "
import pandas as pd
d = pd.read_parquet('results/llm/cluster/$K.parquet')
d.to_csv('results/llm/cluster/$K.csv', index=False)
"
  RESP=max_S_small Rscript code/R/screen_llm.R results/llm/cluster/$K.csv \
      results/llm/cluster_screen_$K 88 1000 > results/llm/cluster_$K.log 2>&1
  rm -f results/llm/cluster/$K.csv
  echo "$(date +%H:%M:%S) $K fait"
done
