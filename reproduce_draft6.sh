#!/bin/sh
# Reproduction pipeline for Draft 6 (article_light_2).
#
# Changes with respect to reproduce_draft5.sh:
#   * the selected tuning of the tail-index screen moves from
#     (a*, b*) = (0.30, 0.15) to (0.35, 0.15);
#   * the Yoshida-Umezu competitor receives a tuning grid search of its
#     own, on the same 200 datasets used to tune our screen, and is run
#     at the best (h, k) FOR EACH MODEL as well as at its paper tuning.
#
# The tuning study of step 1 is unchanged and need not be rerun if
# results/tuning6 is already present: only which cell is reported changes.
#
# Seed streams
#   tuning (both methods)  141000041 + model_index*100019 + r*211
#   comparison datasets    131000021 + cell*100003 + model_index*10007 + r*307
#   tie-breaks             973000019 + cell*100003 + model_index*10007 + r*307
#
# Rough cost on 7 cores: step 2 about 20 min, step 3 several hours.
# Both are checkpointed per block and safe to interrupt and restart.
set -eu

CORES="${CORES:-7}"
NREP_TUNE="${NREP_TUNE:-200}"
NREP_CMP="${NREP_CMP:-1000}"

# 0. Checks.  The YU grid scorer must agree with the reference
#    implementation exactly before it is allowed to choose a tuning.
Rscript tests/test_generate5.R
Rscript tests/test_campaign_equivalence.R
Rscript tests/test_yu_grid.R

# 1. Tuning study for the proposed screen (skip if results/tuning6 exists).
if [ ! -f results/tuning6/summary.csv ]; then
  Rscript code/R/tuning6.R results/tuning6 "$NREP_TUNE" "$CORES"
fi
Rscript code/R/draft5_outputs.R tuning

# 2. Tuning study for the Yoshida-Umezu competitor, same datasets.
Rscript code/R/yu_tuning8.R results/yu_tuning8 "$NREP_TUNE" "$CORES"

# 3. Final comparison: (a*,b*) = (0.35,0.15); YU tuned per model and at
#    its paper tuning.
Rscript code/R/campaign8.R results/campaign8 A1,A2,A3,B1 "$NREP_CMP" "$CORES"

# 4. Tables.
Rscript code/R/draft6_outputs.R comparison

# 5. Real data (aggregated screen; unchanged if the tuning study is unchanged).
Rscript code/R/real_crime5.R results/draft4_real_crime.rds

# 6. Manuscript.
(cd article_light_2 && latexmk -g -pdf -interaction=nonstopmode -halt-on-error main.tex)
