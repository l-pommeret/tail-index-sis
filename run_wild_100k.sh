#!/bin/sh
# Criblage a n = 1e5 : le regime ou 2 n alpha h passe de 71 a 1778.
set -u
cd /people/pommeret/tail-index-sis
PY=/people/pommeret/miniconda3/bin/python
export PATH=/people/pommeret/tail-index-sis/.renv/bin:$PATH
log() { echo "$(date +%H:%M:%S) $*"; }

log "extraction perplexite, n = 100000, GPU 1"
$PY code/py/extract_ppl_same_span.py results/wild/main_docs.parquet \
    results/wild/ppl100k 100000 1 > results/wild/ppl100k.log 2>&1
log "  -> code $?"

log "export en rangs"
$PY code/py/export_ppl_ranks.py results/wild/ppl100k/ppl_cov.parquet \
    results/wild/ppl100k_ranks > results/wild/export100k.log 2>&1
log "  -> code $?"

log "diagnostic de variation de gamma"
$PY code/py/gamma_varies_ppl.py results/wild/ppl100k/ppl_cov.parquet 500 4 0.20 \
    > results/wild/gammavar100k.log 2>&1
log "  -> code $?"

log "criblage"
Rscript code/R/screen_wild.R results/wild/ppl100k_ranks results/wild/screen100k 40 1000 \
    > results/wild/screen100k.log 2>&1
log "  -> code $?"
log "CHAINE 100K TERMINEE"
