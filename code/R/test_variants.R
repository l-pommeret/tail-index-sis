## Do envelope-variation functionals separate better than the envelope level?
## For a few replications of a cell, rank the coordinates by each functional and
## report the realised separation and the worst active rank.  The "mean"
## variant reproduces the current score and is the control: it must match
## score_coordinate_cpp exactly.
## usage: KAPPA=0.20 Rscript code/R/test_variants.R CELL.rds [NREP] [CORES]
suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
sourceCpp("code/src/local_hill_variants.cpp", rebuild = FALSE)
source("code/R/generate3.R")
EPS <- 0.05
args <- commandArgs(trailingOnly = TRUE)
x <- readRDS(args[1L])
nrep <- if (length(args) >= 2L) as.integer(args[2L]) else 8L
cores <- if (length(args) >= 3L) as.integer(args[3L]) else 8L
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))

## direction: TRUE when an active coordinate is expected to give a SMALL value
SMALL_IS_ACTIVE <- c(mean = TRUE, min = TRUE, range = FALSE, sd = FALSE,
                     slope = TRUE)

one <- function(z) {
  n <- z$n; p <- z$p
  d <- simulate_dataset3(n, p, z$rho, z$model, z$seed)
  h <- n^(-BSTAR) / 2; alpha <- n^(-ASTAR)
  V <- vapply(seq_len(p), function(j)
    score_variants_cpp(d$z[, j], d$y, h, alpha, EPS), numeric(6))
  ref <- vapply(seq_len(p), function(j)
    score_coordinate_cpp(d$z[, j], d$y, h, alpha, EPS), numeric(5))["score", ]
  stopifnot(isTRUE(all.equal(unname(V["mean", ]), unname(ref), tolerance = 0)))
  do.call(rbind, lapply(names(SMALL_IS_ACTIVE), function(k) {
    v <- V[k, ]
    s <- if (SMALL_IS_ACTIVE[[k]]) v else -v
    nul <- s[9:p]
    data.frame(variant = k, rep = z$replicate,
               snr = (median(nul) - max(s[1:4])) / sd(nul),
               rmax = max(match(1:4, order(s, seq_len(p)))))
  }))
}
out <- do.call(rbind, parallel::mclapply(x[seq_len(nrep)], one, mc.cores = cores))
z1 <- x[[1]]
thr <- qnorm(1 - 1 / (z1$p - 8))
cat(sprintf("\n%s n=%d p=%d rho=%.2f, %d replications, a*=%.2f b*=%.2f\n",
            z1$model, z1$n, z1$p, z1$rho, nrep, ASTAR, BSTAR))
cat(sprintf("seuil de competition : minimum de %d scores nuls a ~%.2f sigma\n\n",
            z1$p - 8L, thr))
cat(sprintf("%-8s %8s %8s %8s %8s %8s\n",
            "variante", "SNR moy", "SNR med", "Rmax med", "Rmax max", "P(<=20)"))
for (k in names(SMALL_IS_ACTIVE)) {
  w <- out[out$variant == k, ]
  cat(sprintf("%-8s %8.2f %8.2f %8.0f %8.0f %8.2f\n", k, mean(w$snr),
              median(w$snr), median(w$rmax), max(w$rmax), mean(w$rmax <= 20)))
}
