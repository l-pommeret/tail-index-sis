#!/bin/sh
set -u
cd /people/pommeret/tail-index-sis
export PATH=/people/pommeret/tail-index-sis/.renv/bin:$PATH
# attendre la fin du balayage a alpha fixe pour ne pas se disputer les coeurs
while pgrep -f "wild_smalln.R results/wild/ppl100k_ranks results/wild/nsweep_all/nsweep_all.rds results/wild/smalln " > /dev/null 2>&1; do sleep 20; done
echo "$(date +%H:%M:%S) petits n, convention a = 0.35"
Rscript code/R/wild_smalln.R results/wild/ppl100k_ranks \
  results/wild/nsweep_all/nsweep_all.rds results/wild/smalln_a035 40 20 fixed_a \
  > results/wild/smalln_a035.log 2>&1
echo "$(date +%H:%M:%S)   -> code $?"
