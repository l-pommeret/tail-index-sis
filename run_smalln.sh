#!/bin/sh
set -u
cd /people/pommeret/tail-index-sis
export PATH=/people/pommeret/tail-index-sis/.renv/bin:$PATH
echo "$(date +%H:%M:%S) petits n, 7 regles, 20 tirages"
Rscript code/R/wild_smalln.R results/wild/ppl100k_ranks \
  results/wild/nsweep_all/nsweep_all.rds results/wild/smalln 40 20 \
  > results/wild/smalln.log 2>&1
echo "$(date +%H:%M:%S)   -> code $?"
