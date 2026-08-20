#!/bin/sh
set -u
cd /people/pommeret/tail-index-sis
export PATH=/people/pommeret/tail-index-sis/.renv/bin:$PATH
echo "$(date +%H:%M:%S) petits n, convention a = 0.35"
Rscript code/R/wild_smalln.R results/wild/ppl100k_ranks \
  results/wild/nsweep_all/nsweep_all.rds results/wild/smalln_a035 40 20 fixed_a \
  > results/wild/smalln_a035.log 2>&1
echo "$(date +%H:%M:%S) A035 FINI code $?"
