#!/bin/sh
set -u
cd /people/pommeret/tail-index-sis
export PATH=/people/pommeret/tail-index-sis/.renv/bin:$PATH
echo "$(date +%H:%M:%S) grille (n, a)"
Rscript code/R/wild_grid_na.R results/wild/ppl100k_ranks \
  results/wild/nsweep_all/nsweep_all.rds results/wild/grid_na 40 \
  > results/wild/grid_na.log 2>&1
echo "$(date +%H:%M:%S) GRID FINI code $?"
