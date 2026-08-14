#!/bin/sh
# Reproduction pipeline for Draft 2.
# Requires: R with Rcpp and quantreg, a C++ toolchain, LaTeX (latexmk).
set -eu

# 1. Tuning heatmaps: 4 models x 42 (a,b) cells x 200 replications.
Rscript code/R/run_draft2b.R tuning results/draft2/tuning_cells 7   # checkpointed, resumable

# 2. Fix (a*, b*) from the common robust region, then one-factor
#    sensitivity on Model 1 (independent seed stream).
# Selected from the common robust plateau of the four heatmaps:
ASTAR=${ASTAR:-0.30}
BSTAR=${BSTAR:-0.10}
export ASTAR BSTAR
Rscript code/R/run_draft2b.R sensitivity results/draft2/sensitivity_cells 7

# 3. Three-method comparison on common datasets (independent seed stream).
Rscript code/R/run_draft2b.R comparison results/draft2/comparison_cells 7

# 3b. Competitor-fairness supplement: Yoshida-Umezu bandwidth grid.
Rscript code/R/run_yu_grid2.R results/draft2/yu_grid_cells 7

# 4. Real data: Communities and Crime (Unnormalized), UCI id 211.
#    data/crime/CommViolPredUnnormalizedData.txt from
#    https://archive.ics.uci.edu/static/public/211/communities+and+crime+unnormalized.zip
Rscript code/R/real_crime.R "$ASTAR" "$BSTAR" results/draft2/real_crime.rds
Rscript code/R/real_crime_outputs.R

# 5. Manuscript.
(cd manuscript && latexmk -g -pdf -interaction=nonstopmode -halt-on-error main.tex)
shasum -a 256 results/draft2/*.rds results/draft2/*.csv > results/draft2/SHA256SUMS.txt
