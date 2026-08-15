#!/bin/sh
set -u
cd /people/pommeret/tail-index-sis
PY=/people/pommeret/miniconda3/bin/python
export PATH=/people/pommeret/tail-index-sis/.renv/bin:$PATH
for C in c03 c06 c18 c24 c3 acp; do
  K=1024
  $PY code/py/export_ranks.py results/llm/internal.parquet results/llm/ranks $C $K \
      >> results/llm/export.log 2>&1
  Rscript code/R/screen_internal.R results/llm/ranks $C results/llm/screen_internal 88 1000 FALSE \
      > results/llm/screen_$C.log 2>&1
  echo "$(date +%H:%M:%S) $C fait"
  rm -f results/llm/ranks/${C}_ranks.bin
done
