## Does the (a, b) tuning chosen at n=2000, p=1000 still separate at n=5000,
## p=2000?  For a handful of replications of one cell, sweep (a, b) and report
## the realised separation in units of the null spread -- the quantity that has
## to beat the minimum of p-4 null scores -- together with the resulting ranks.
## usage: KAPPA=0.20 Rscript code/R/diag_tuning_snr.R CELL.rds [NREP] [CORES]
suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R")
EPS <- 0.05
args <- commandArgs(trailingOnly = TRUE)
x <- readRDS(args[1L])
nrep <- if (length(args) >= 2L) as.integer(args[2L]) else 3L
cores <- if (length(args) >= 3L) as.integer(args[3L]) else 3L

AS <- c(0.30, 0.40, 0.50)
BS <- c(0.10, 0.25, 0.40)
reps <- x[seq_len(nrep)]

one <- function(z) {
  n <- z$n; p <- z$p
  d <- simulate_dataset3(n, p, z$rho, z$model, z$seed)
  do.call(rbind, lapply(AS, function(a) do.call(rbind, lapply(BS, function(b) {
    h <- n^(-b) / 2; alpha <- n^(-a)
    S <- vapply(seq_len(p), function(j)
      score_coordinate_cpp(d$z[, j], d$y, h, alpha, EPS), numeric(5))
    sc <- S["score", ]
    nul <- sc[9:p]
    data.frame(a = a, b = b, rep = z$replicate,
      k_over_m = mean(S["mean_local_k", ]) / mean(S["mean_local_n", ]),
      indep_windows = (1 - 2 * EPS) / (2 * h),
      delta = median(nul) - max(sc[1:4]),
      sigma = sd(nul),
      snr = (median(nul) - max(sc[1:4])) / sd(nul),
      rmax = max(match(1:4, order(sc, seq_len(p)))),
      nonfinite = sum(!is.finite(sc)))
  }))))
}
out <- do.call(rbind, parallel::mclapply(reps, one, mc.cores = cores))
ag <- aggregate(cbind(k_over_m, indep_windows, snr, rmax) ~ a + b, out, mean)
ag <- ag[order(ag$a, ag$b), ]
z1 <- x[[1]]
cat(sprintf("\n%s n=%d p=%d rho=%.2f, %d replications ; reference du papier a*=0.30 b*=0.10\n",
            z1$model, z1$n, z1$p, z1$rho, nrep))
cat(sprintf("seuil de competition : le minimum de %d scores nuls est a environ %.2f sigma sous leur mediane\n\n",
            z1$p - 8L, qnorm(1 - 1 / (z1$p - 8))))
cat(sprintf("%5s %5s %8s %9s %8s %8s\n", "a", "b", "k/m", "fenetres", "SNR", "Rmax"))
for (i in seq_len(nrow(ag)))
  cat(sprintf("%5.2f %5.2f %7.1f%% %9.1f %8.2f %8.1f\n", ag$a[i], ag$b[i],
              100 * ag$k_over_m[i], ag$indep_windows[i], ag$snr[i], ag$rmax[i]))
