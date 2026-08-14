#!/bin/sh
# LEGACY: Draft-1 pipeline. Kept for the archived Draft-1 results only;
# the Draft-2 manuscript is reproduced by reproduce_draft2.sh.
set -eu

for name in correctness sensitivity scaling main_dimension dependence models extension_p10000
do
  Rscript code/R/audit_experiment.R "results/${name}.rds"
  Rscript code/R/summarize_experiment.R "results/${name}.rds" "results/tables/${name}.csv"
  Rscript code/R/extended_metrics.R "results/${name}.rds" "results/tables/${name}_by_d.csv"
done
Rscript code/R/audit_experiment.R results/alpha_large_p5000.rds
Rscript code/R/summarize_experiment.R results/alpha_large_p5000.rds results/tables/alpha_large_p5000.csv
Rscript code/R/extended_metrics.R results/alpha_large_p5000.rds results/tables/alpha_large_p5000_by_d.csv
Rscript code/R/audit_comparators.R results/comparators.rds
Rscript code/R/summarize_comparators.R results/comparators.rds results/tables/comparators.csv
Rscript code/R/audit_yu_tuning_grid.R results/yu_tuning_grid.rds
Rscript code/R/audit_real_online_news.R results/real_online_news/online_news_application.rds
Rscript code/R/final_experiment_outputs.R results results/final
Rscript code/R/canonical_cell.R results results/final
Rscript code/R/make_final_figures.R results/final
Rscript code/R/make_population_figure.R manuscript/figures/population_profiles.pdf
Rscript code/R/make_latex_tables.R results/tables manuscript/tables
cp results/final/sure_screening_vs_p.pdf manuscript/figures/main_dimension.pdf
cp results/final/sensitivity_heatmap.pdf manuscript/figures/sensitivity.pdf
cp results/real_online_news/online_news_application.pdf manuscript/figures/online_news_application.pdf
Rscript tests/test_core.R
shasum -a 256 results/*.rds results/yu_tuning_grid.csv results/tables/*.csv results/final/* > results/SHA256SUMS.txt
(cd manuscript && latexmk -g -pdf -interaction=nonstopmode -halt-on-error main.tex)
