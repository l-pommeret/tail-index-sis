## Replay the worst replication of a cell and dissect the score landscape:
## where do the active coordinates sit, what does the null look like, and are
## the local windows healthy?  Costs one dataset, so it runs on a single core.
## usage: KAPPA=0.20 Rscript code/R/diag_failure.R CELL.rds
suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R")
EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))

f <- commandArgs(trailingOnly = TRUE)[1L]
x <- readRDS(f)
rmax <- vapply(x, function(z) max(z$ranks_ours), numeric(1))
for (tag in c("pire", "median")) {
  i <- if (tag == "pire") which.max(rmax) else which.min(abs(rmax - median(rmax)))
  z <- x[[i]]
  n <- z$n; p <- z$p
  d <- simulate_dataset3(n, p, z$rho, z$model, z$seed)
  h <- n^(-BSTAR) / 2; alpha <- n^(-ASTAR)
  S <- vapply(seq_len(p), function(j)
    score_coordinate_cpp(d$z[, j], d$y, h, alpha, EPS), numeric(5))
  sc <- S["score", ]
  cat(sprintf("\n=== %s : replication %d, seed %d, Rmax = %d  (%s n=%d p=%d rho=%.2f)\n",
              tag, z$replicate, z$seed, max(z$ranks_ours), z$model, n, p, z$rho))
  cat(sprintf("  gamma vrai : min %.3f  median %.3f  max %.3f\n",
              min(d$gamma), median(d$gamma), max(d$gamma)))
  cat(sprintf("  fenetre locale : h = %.4f (largeur %.1f%% des donnees), m = %.0f, k = %.0f, fenetres independantes ~ %.1f\n",
              h, 200 * h, mean(S["mean_local_n", ]), mean(S["mean_local_k", ]),
              (1 - 2 * EPS) / (2 * h)))
  cat(sprintf("  scores non finis : %d ; taux de sous-peuplement moyen : %.3f (max %.3f)\n",
              sum(!is.finite(sc)), mean(S["under_rate", ]), max(S["under_rate", ])))
  act <- sc[1:4]; sca <- sc[5:8]; nul <- sc[9:p]
  cat(sprintf("  actives  X1..X4 : %s\n", paste(sprintf("%.4f", act), collapse = " ")))
  cat(sprintf("  echelle  X5..X8 : %s\n", paste(sprintf("%.4f", sca), collapse = " ")))
  qn <- quantile(nul, c(0, .001, .01, .05, .5, .95, 1))
  cat(sprintf("  nulles (n=%d)   : min %.4f  q0.1%% %.4f  q1%% %.4f  q5%% %.4f  med %.4f  max %.4f\n",
              length(nul), qn[1], qn[2], qn[3], qn[4], qn[5], qn[7]))
  cat(sprintf("  ecart-type des nulles %.4f ; separation (med nulles - pire active) = %.4f = %.1f sigma\n",
              sd(nul), qn[5] - max(act), (qn[5] - max(act)) / sd(nul)))
  cat(sprintf("  rangs actifs : %s\n", paste(match(1:4, order(sc, seq_len(p))),
                                             collapse = " ")))
}
