## Draft-3 M2 calibration pilot (model calibration only; results never enter
## the final comparison).  Pre-specified protocol:
##   kappa in {0.15, 0.25, 0.35, 0.50, 0.70}, 80 replications each,
##   baseline design (n, p, rho) = (2000, 1000, 0.25),
##   common tuning (a*, b*) = (0.30, 0.10),
##   dedicated pilot seed stream 7000003 + k_index*100057 + r*503,
##   disjoint from the tuning (2000003+...), sensitivity (3000017+...),
##   comparison (4000037+...) and real-data streams.
## Selection rule, fixed before inspection: choose the SMALLEST kappa with
##   (i) Sure-20(Tail) >= 0.5,
##   (ii) E|A4_hat(Q95) ∩ A_scale| clearly above E|A4_hat(Tail) ∩ A_scale|,
##   (iii) E|A4_hat(Tail) ∩ A_gamma| > E|A4_hat(Q95) ∩ A_gamma|,
##   (iv) visible scale effect: E|A4_hat(Q95) ∩ A_scale| >= 1.
## usage: KAPPA ignored here; Rscript code/R/pilot_m2.R OUT.rds CORES
args <- commandArgs(trailingOnly = TRUE)
out_path <- args[1L]; cores <- as.integer(args[2L])
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R"); source("code/R/qa_sis.R")

## One pre-authorized grid extension: the initial grid bracketed the
## target between 0.15 (too weak a scale effect) and 0.25 (tail signal
## degraded); kappa = 0.20 is the single midpoint added.
KGRID <- c(0.15, 0.25, 0.35, 0.50, 0.70, 0.20)
n <- 2000L; p <- 1000L; rho <- 0.25
h <- n^(-0.10) / 2; alpha <- n^(-0.30)

res <- list()
for (ki in seq_along(KGRID)) {
  kap <- KGRID[ki]
  assign("KAPPA", kap, envir = globalenv())
  out <- parallel::mclapply(1:80, function(r) {
    seed <- 7000003L + ki * 100057L + r * 503L
    d <- simulate_dataset3(n, p, rho, "M2", seed)
    sc <- vapply(seq_len(p), function(j)
      score_coordinate_cpp(d$z[, j], d$y, h, alpha, 0.05), numeric(5))["score", ]
    o_t <- order(sc, seq_len(p), na.last = TRUE)
    q95 <- qa_sis_scores(
      apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (n + 1)),
      d$y, tau = 0.95)
    o_q <- order(-q95, seq_len(p))
    list(kappa = kap, replicate = r, seed = seed,
         rt = match(1:4, o_t), rq = match(1:4, o_q),
         top4_t = o_t[1:4], top4_q = o_q[1:4])
  }, mc.cores = cores, mc.preschedule = FALSE)
  res[[as.character(kap)]] <- out
  saveRDS(res, out_path, compress = "xz")  # checkpoint after each kappa
  cat(format(Sys.time(), "%H:%M:%S"), "pilot kappa", kap, "done\n")
}
saveRDS(res, out_path, compress = "xz")

for (kap in names(res)) {
  x <- res[[kap]]
  rt <- sapply(x, function(z) max(z$rt)); rq <- sapply(x, function(z) max(z$rq))
  cat(sprintf(
    "kappa=%s Tail: S4=%.2f S20=%.2f Rmax=%.0f | Q95: S4=%.2f S20=%.2f Rmax=%.0f | top4 Tail g=%.2f s=%.2f | Q95 g=%.2f s=%.2f\n",
    kap, mean(rt <= 4), mean(rt <= 20), mean(rt),
    mean(rq <= 4), mean(rq <= 20), mean(rq),
    mean(sapply(x, function(z) sum(z$top4_t %in% 1:4))),
    mean(sapply(x, function(z) sum(z$top4_t %in% 5:8))),
    mean(sapply(x, function(z) sum(z$top4_q %in% 1:4))),
    mean(sapply(x, function(z) sum(z$top4_q %in% 5:8)))))
}
