#!/bin/sh
set -u
cd /people/pommeret/tail-index-sis
export PATH=/people/pommeret/tail-index-sis/.renv/bin:$PATH
echo "$(date +%H:%M:%S) balayage en n"
Rscript code/R/wild_nsweep.R results/wild/ppl100k_ranks results/wild/nsweep 40 400 \
  > results/wild/nsweep.log 2>&1
echo "$(date +%H:%M:%S)   -> code $?"
