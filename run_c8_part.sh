#!/bin/sh
# Une part de campaign8. Les blocs deja presents sont sautes (checkpoints),
# donc les deux machines peuvent travailler sur le meme repertoire tant que
# leurs listes de modeles sont disjointes.
set -u
cd /people/pommeret/tail-index-sis
export PATH=/people/pommeret/tail-index-sis/.renv/bin:$PATH
MODELS="$1"; C="$2"; TAG="$3"
echo "$(date +%H:%M:%S) [$TAG] campaign8 $MODELS sur $C coeurs"
C8_BLOCK=50 Rscript code/R/campaign8.R results/campaign8 "$MODELS" 1000 "$C" \
    > "results/campaign8_$TAG.log" 2>&1
echo "$(date +%H:%M:%S) [$TAG] fini, code $?"
