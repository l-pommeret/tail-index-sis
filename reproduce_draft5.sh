#!/bin/sh
# Reproduction pipeline for Draft 5.
# Models A1/A2/A3/B1 are defined in code/R/generate5.R; the unit checks in
# tests/test_generate5.R verify that the full and streaming generators
# implement the same distribution (including the corrected continuation of
# the AR(1) chain past the B1 proxy block), and
# tests/test_campaign_equivalence.R verifies that A1 and A2 reproduce the
# campaign5 M1 and M3 datasets exactly, which justifies relabelling those
# stored results instead of rerunning them.
# Seed streams: tuning6 141000041+m*100019+r*211 (common random numbers
# across all 42 cells); comparison 131000021+cell*100003+model*10007+r*307
# (campaign5 layout); real data as in code/R/real_crime5.R.
set -eu

Rscript tests/test_generate5.R
Rscript tests/test_campaign_equivalence.R

# 1. Tuning study: 200 common-random-number replications per model,
#    every (a,b) cell evaluated on the same datasets.
Rscript code/R/tuning6.R results/tuning6 200 7
Rscript code/R/draft5_outputs.R tuning

# 2. Final comparison: A3 and B1 rerun under generate5.R; A1 and A2 reused
#    from results/campaign5 (equivalence verified above).
Rscript code/R/campaign6.R results/campaign6 A3,B1 1000 7
Rscript code/R/draft5_outputs.R comparison

# 3. Real data (aggregated screen; unchanged if tuning is unchanged).
Rscript code/R/real_crime5.R results/draft4_real_crime.rds

# 4. Manuscript.
(cd manuscript && latexmk -g -pdf -interaction=nonstopmode -halt-on-error main.tex)
