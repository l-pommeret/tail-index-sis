#!/bin/sh
set -u
cd /people/pommeret/tail-index-sis
export PATH=/people/pommeret/tail-index-sis/.renv/bin:$PATH
echo "$(date +%H:%M:%S) criblage CORRIGE, a=0.20 (alpha=0.10), n=1e5"
WILD_A=0.20 WILD_B=0.15 Rscript code/R/screen_wild.R results/wild/ppl100k_ranks \
    results/wild/screen100k_fix 40 1000 > results/wild/screen100k_fix.log 2>&1
echo "$(date +%H:%M:%S)   -> code $?"
echo "$(date +%H:%M:%S) criblage CORRIGE, a=0.35 (reglage manuscrit)"
WILD_A=0.35 WILD_B=0.15 Rscript code/R/screen_wild.R results/wild/ppl100k_ranks \
    results/wild/screen100k_fix_a035 40 1000 > results/wild/screen100k_fix_a035.log 2>&1
echo "$(date +%H:%M:%S)   -> code $?"
echo "$(date +%H:%M:%S) RESCREEN TERMINE"
