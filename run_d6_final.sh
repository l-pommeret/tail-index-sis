#!/bin/sh
# Attend la fin des deux parts de campaign8, puis enchaine les etapes 4 et 5.
set -u
cd /people/pommeret/tail-index-sis
export PATH=/people/pommeret/tail-index-sis/.renv/bin:$PATH
done_p() { grep -q "fini, code" "$1" 2>/dev/null; }
while true; do
  if done_p results/part_giles.log && done_p results/part_baudelaire.log; then break; fi
  sleep 30
done
echo "$(date +%H:%M:%S) les deux parts sont finies ; $(ls results/campaign8/*.rds | wc -l) blocs"
echo "$(date +%H:%M:%S) resume complet de campaign8"
C8_BLOCK=50 Rscript code/R/campaign8.R results/campaign8 A1,A2,A3,B1 1000 8 > results/campaign8_resume.log 2>&1
echo "  -> code $?"
echo "$(date +%H:%M:%S) etape 4 : tables de comparaison"
Rscript code/R/draft6_outputs.R comparison > results/d6_cmp_tables.log 2>&1
echo "  -> code $?"
echo "$(date +%H:%M:%S) etape 5 : donnees reelles"
Rscript code/R/real_crime5.R results/draft4_real_crime.rds > results/d6_real_crime.log 2>&1
echo "  -> code $?"
echo "$(date +%H:%M:%S) DRAFT 6 TOUT TERMINE"
