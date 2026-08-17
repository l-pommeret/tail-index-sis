#!/bin/sh
set -u
cd /people/pommeret/tail-index-sis
export PATH=/people/pommeret/tail-index-sis/.renv/bin:$PATH
echo "$(date +%H:%M:%S) application wild : les sept regles"
Rscript code/R/real_wild.R results/wild/ppl100k_ranks \
  results/wild/ppl100k/ppl_cov_gammavar.csv results/wild/real_wild.rds 40 \
  > results/wild/real_wild.log 2>&1
echo "$(date +%H:%M:%S)   -> code $?"
echo "$(date +%H:%M:%S) REAL WILD TERMINE"
