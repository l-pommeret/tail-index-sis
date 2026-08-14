#!/bin/sh
# Reproduction pipeline for Draft 3.
# Requires: R with Rcpp and quantreg, a C++ toolchain, LaTeX (latexmk).
# Provenance: pilot stream 7000003+k*100057+r*503 (calibration only);
# tuning-M2 stream 2200003+i*100003+r*211; comparison stream
# 4000037+i*100043+r*211 (M1/M3/M4 datasets identical to Draft 2);
# real-data half-samples seeds 500001..500200, tie audit 900001..900020.
set -eu

export KAPPA=${KAPPA:-0.20}   # fixed by the pre-specified pilot (below)
export ASTAR=${ASTAR:-0.30}
export BSTAR=${BSTAR:-0.10}

# 0. M2 calibration pilot (model calibration only; never enters final results).
Rscript code/R/pilot_m2.R results/draft3/pilot_m2.rds 7
Rscript code/R/pilot_table.R

# 1. Tuning: M1/M3/M4 heatmaps are unchanged from Draft 2
#    (results/draft2/tuning_cells); recompute the M2 heatmap only.
Rscript code/R/run_draft3.R tuning_m2 results/draft3/tuning_m2_cells 7
Rscript code/R/draft3_outputs.R tuning

# 2. One-factor sensitivity (Model M1, unchanged from Draft 2):
#    results/draft2/sensitivity_cells.

# 3. Method comparison, all four models, quantile SIS at
#    tau in {.90,.95,.975,.99}, on common datasets; YU bandwidth grid for M2.
Rscript code/R/run_draft3.R comparison results/draft3/comparison_cells 7
Rscript code/R/run_draft3.R yugrid_m2 results/draft3/yugrid_m2_cells 7
Rscript code/R/draft3_outputs.R comparison

# 4. Real data (harmonized median-imputed preprocessing, Hill + Pareto QQ).
Rscript code/R/real_crime.R "$ASTAR" "$BSTAR" results/draft3/real_crime.rds
Rscript code/R/draft3_outputs.R real

# 5. Manuscript.
(cd manuscript && latexmk -g -pdf -interaction=nonstopmode -halt-on-error main.tex)
shasum -a 256 results/draft3/*.csv results/draft3/*.rds > results/draft3/SHA256SUMS.txt
