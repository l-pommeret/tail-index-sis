## Is the bimodality of Rmax a property of the DATASET or of the ESTIMATOR?
##
## The distinction decides what can fix it.  If a replication fails because the
## estimator was unlucky on that draw, averaging over tunings or over subsamples
## of the SAME data will rescue it.  If the realised separation collapses for
## that dataset, every tuning and every subsample fails together and no internal
## aggregation can help.
##
## For each replication this computes:
##   - the baseline rank at (a*, b*);
##   - ranks over a 3x3 tuning grid, aggregated by median rank;
##   - ranks over B half-sample refits, aggregated by median rank;
## and records, per replication, the share of tunings and of half-samples that
## fail, which is the quantity that separates the two explanations: a share near
## 0 or 1 means the failure is carried by the dataset, an intermediate share
## means it is estimator noise that aggregation can average away.
##
## usage: KAPPA=0.20 Rscript code/R/test_stability.R CORES CELL.rds [...]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R")

args <- commandArgs(trailingOnly = TRUE)
CORES <- as.integer(args[1L])
cells <- args[-1L]
dir.create("results/stability", showWarnings = FALSE, recursive = TRUE)
EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
AGRID <- c(0.30, 0.35, 0.40)
BGRID <- c(0.05, 0.10, 0.15)
NSUB <- 10L                      # half-sample refits per replication

score_vec <- function(z, y, a, b) {
  n <- nrow(z); p <- ncol(z)
  h <- n^(-b) / 2; alpha <- n^(-a)
  vapply(seq_len(p), function(j)
    score_coordinate_cpp(z[, j], y, h, alpha, EPS), numeric(5))["score", ]
}
rank_of <- function(s, p) rank(s, ties.method = "first", na.last = TRUE)

one <- function(z0) {
  n <- z0$n; p <- z0$p
  d <- simulate_dataset3(n, p, z0$rho, z0$model, z0$seed)
  base <- score_vec(d$z, d$y, ASTAR, BSTAR)
  rmax_base <- max(match(1:4, order(base, seq_len(p), na.last = TRUE)))

  ## --- tuning grid -----------------------------------------------------------
  R <- matrix(NA_real_, p, length(AGRID) * length(BGRID))
  rmax_tune <- numeric(ncol(R)); ii <- 0L
  for (a in AGRID) for (b in BGRID) {
    ii <- ii + 1L
    s <- score_vec(d$z, d$y, a, b)
    R[, ii] <- rank_of(s, p)
    rmax_tune[ii] <- max(match(1:4, order(s, seq_len(p), na.last = TRUE)))
  }
  agg_t <- apply(R, 1L, median)
  rmax_agg_tune <- max(match(1:4, order(agg_t, seq_len(p))))

  ## --- half-sample refits ----------------------------------------------------
  set.seed(z0$seed + 77L)
  S <- matrix(NA_real_, p, NSUB)
  rmax_sub <- numeric(NSUB)
  for (bb in seq_len(NSUB)) {
    idx <- sample.int(n, floor(n / 2))
    s <- score_vec(d$z[idx, , drop = FALSE], d$y[idx], ASTAR, BSTAR)
    S[, bb] <- rank_of(s, p)
    rmax_sub[bb] <- max(match(1:4, order(s, seq_len(p), na.last = TRUE)))
  }
  agg_s <- apply(S, 1L, median)
  rmax_agg_sub <- max(match(1:4, order(agg_s, seq_len(p))))

  list(model = z0$model, n = n, p = p, rho = z0$rho,
       rmax_base = rmax_base, rmax_agg_tune = rmax_agg_tune,
       rmax_agg_sub = rmax_agg_sub,
       frac_tune_fail = mean(rmax_tune > 20), frac_sub_fail = mean(rmax_sub > 20),
       med_sub = median(rmax_sub))
}

res <- list()
for (f in cells) {
  x <- readRDS(f)
  out <- parallel::mclapply(x, one, mc.cores = CORES, mc.preschedule = FALSE)
  bad <- vapply(out, function(z) inherits(z, "try-error") || is.null(z),
                logical(1))
  if (any(bad)) stop(f, ": ", sum(bad), " failed: ", as.character(out[bad][[1]]))
  res[[f]] <- out
  cat(format(Sys.time(), "%H:%M:%S"), basename(f), "done\n"); flush.console()
}
saveRDS(res, "results/stability/raw.rds", compress = "xz")

fm <- function(v, k = 3) formatC(v, format = "f", digits = k)
cat("\nStabilite de Rmax : agregation sur reglages et sur demi-echantillons\n\n")
for (f in cells) {
  out <- res[[f]]; z1 <- out[[1]]
  g <- function(w) vapply(out, `[[`, numeric(1), w)
  cat(sprintf("=== %s  n=%d p=%d rho=%.2f, %d replications\n",
              z1$model, z1$n, z1$p, z1$rho, length(out)))
  cat(sprintf("%-28s %8s %8s %8s\n", "regle", "Sure-4", "Sure-10", "Sure-20"))
  for (nm in c("rmax_base", "rmax_agg_tune", "rmax_agg_sub")) {
    v <- g(nm)
    cat(sprintf("%-28s %8s %8s %8s\n",
                switch(nm, rmax_base = "score seul (a*, b*)",
                       rmax_agg_tune = "rang median sur 9 reglages",
                       rmax_agg_sub = "rang median sur 10 demi-ech."),
                fm(mean(v <= 4)), fm(mean(v <= 10)), fm(mean(v <= 20))))
  }
  ## do the failures coincide?  split replications by whether the baseline failed
  ok <- g("rmax_base") <= 20
  cat(sprintf("\n  part des 9 reglages qui echouent   : %.2f quand le score seul reussit, %.2f quand il echoue\n",
              mean(g("frac_tune_fail")[ok]), mean(g("frac_tune_fail")[!ok])))
  cat(sprintf("  part des 10 demi-ech. qui echouent : %.2f quand le score seul reussit, %.2f quand il echoue\n\n",
              mean(g("frac_sub_fail")[ok]), mean(g("frac_sub_fail")[!ok])))
}
