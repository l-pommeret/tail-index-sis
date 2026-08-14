# Tail-index sure independence screening

Research companion for **Tail-Index Sure Independence Screening for
Ultra-High-Dimensional Heavy-Tailed Regression** (Draft 3).

The paper develops a coordinatewise projected-tail-index screen: after
marginal probability-integral transformation, conditioning on one coordinate
turns the conditional tail-index surface into its upper envelope along the
corresponding fibre. Inactive coordinates give flat envelopes at the global
maximum; detectable active coordinates depress them. The screen ranks
coordinates by empirical-rank local Hill averages of these envelopes, with a
finite-sample uniform bound over a diverging number of coordinates, sure
screening under `log(pn) = o(n α h Δ²_min)`, an `O(p n log n)` implementation,
a four-model simulation study built on the slowly varying factors of Gardes
and Podgorny, a quantile-level sensitivity study, and an application to
violent-crime rates across 1,993 U.S. communities (p = 99 ≈ effective tail
sample).

The Lean companion verifies the deterministic bridge from uniform score error
to separation, sure screening and exact recovery:
<https://github.com/l-pommeret/tail-index-sis-lean>.

## Layout

- `manuscript/` — LaTeX sources, figures, tables, `main.pdf`.
- `code/R/`, `code/src/` — estimator (Rcpp Fenwick-tree local Hill),
  simulation drivers, competitor implementations (Yoshida–Umezu
  conditional-Pickands screen; He–Wang–Hong quantile-adaptive SIS),
  M2 calibration pilot, real-data analysis, table/figure generation.
- `results/draft2/`, `results/draft3/` — replicate-level checkpoints
  (`*_cells/`), summary CSVs, and `SHA256SUMS.txt`. Every number quoted in
  the manuscript traces to these files.
- `data/crime/` — Communities and Crime (Unnormalized), UCI Machine Learning
  Repository id 211 (CC BY 4.0; Redmond & Baveja 2002). Included verbatim for
  reproducibility; `reproduce_draft3.sh` documents the download URL.
- `reproduce_draft3.sh` — full pipeline (pilot → tuning → comparison →
  application → manuscript), with seed-stream provenance in the header.
  `reproduce_draft2.sh` reproduces the Draft-2 experiments that Draft 3
  carries over unchanged (M1/M3/M4 tuning heatmaps, one-factor sensitivity).
- `REVISION_LOG.md` — changes across drafts, including the pre-registered
  M2 calibration protocol and the review-driven corrections.

## Requirements

R (≥ 4.5) with `Rcpp` and `quantreg`, a C++ toolchain, TeX Live with
`latexmk`. Run all commands from the repository root.

## Reproduction

```sh
./reproduce_draft3.sh
```

Individual stages can be rerun independently; all drivers are checkpointed
per cell and resume automatically. Monte Carlo seed streams for the pilot,
tuning, sensitivity, comparison and application are disjoint and documented
in `reproduce_draft3.sh` and `code/R/pilot_m2.R`.
