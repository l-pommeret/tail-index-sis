#!/bin/sh
# Criblage du mode libre des que l'extraction est ecrite : c'est le test M2,
# ou la longueur varie et doit se comporter en nuisance d'echelle pure.
set -u
cd /people/pommeret/tail-index-sis
export PATH=/people/pommeret/tail-index-sis/.renv/bin:$PATH
until [ -f results/llm/main_libre.parquet ]; do sleep 60; done
sleep 10
/people/pommeret/miniconda3/bin/python -c "
import pandas as pd
d = pd.read_parquet('results/llm/main_libre.parquet')
d.to_csv('results/llm/main_libre.csv', index=False)
print('CSV', d.shape)
" > results/llm/libre_csv.log 2>&1
echo "$(date +%H:%M:%S) csv ecrit"
Rscript code/R/screen_llm.R results/llm/main_libre.csv results/llm/screen_libre 88 1000 \
    > results/llm/screen_libre.log 2>&1
echo "$(date +%H:%M:%S) criblage libre termine"
