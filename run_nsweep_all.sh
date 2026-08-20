#!/bin/sh
set -u
cd /people/pommeret/tail-index-sis
export PATH=/people/pommeret/tail-index-sis/.renv/bin:$PATH
echo "$(date +%H:%M:%S) balayage en n, les sept regles"
Rscript code/R/wild_nsweep_all.R results/wild/ppl100k_ranks results/wild/nsweep_all 40 \
  > results/wild/nsweep_all.log 2>&1
echo "$(date +%H:%M:%S)   -> code $?"
