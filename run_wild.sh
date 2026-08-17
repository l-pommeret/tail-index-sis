#!/bin/sh
# Etage 1 : corpus principal FineWeb a n = 1e5 et reponse Y = exp(max NLL).
set -u
cd /people/pommeret/tail-index-sis
PY=/people/pommeret/miniconda3/bin/python
log() { echo "$(date +%H:%M:%S) $*"; }
log "tirage 120000 documents"
$PY code/py/sample_fineweb.py results/wild/main_docs.parquet 120000 > results/wild/sample.log 2>&1
log "  -> code $?"
log "extraction NLL sur GPU 1"
$PY code/py/extract_nll.py results/wild/main_docs.parquet results/wild/main_nll.parquet \
    100000 EleutherAI/pythia-410m 1 > results/wild/nll.log 2>&1
log "  -> code $?"
log "ETAGE 1 TERMINE"
