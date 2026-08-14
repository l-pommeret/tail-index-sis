## Does second-order bias reduction restore the separation the local Hill loses?
## Reports, per rule, the screening performance and the realised separation in
## units of the null spread -- the quantity that must beat the minimum of p-4
## null scores.  The plain Hill on the SAME evaluation grid is the control, so
## any difference is the bias correction and not the grid.
## usage: KAPPA=0.20 Rscript code/R/test_bias_reduced.R CORES CELL.rds [...]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
sourceCpp("code/src/bias_reduced_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R")

args <- commandArgs(trailingOnly = TRUE)
CORES <- as.integer(args[1L])
cells <- args[-1L]
dir.create("results/bias_reduced", showWarnings = FALSE, recursive = TRUE)
EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
RHOS <- c(-1.0, -2.0, -3.0, -4.0)
RULES <- c("Hill (score campagne)", "Hill (grille 25)",
           paste0("biais reduit rho=", RHOS))

one <- function(z) {
  n <- z$n; p <- z$p
  d <- simulate_dataset3(n, p, z$rho, z$model, z$seed)
  h <- n^(-BSTAR) / 2; alpha <- n^(-ASTAR)
  ref <- vapply(seq_len(p), function(j)
    score_coordinate_cpp(d$z[, j], d$y, h, alpha, EPS), numeric(5))["score", ]
  sc <- list(ref)
  for (i in seq_along(RHOS)) {
    V <- vapply(seq_len(p), function(j)
      score_bias_reduced_cpp(d$z[, j], d$y, h, alpha, RHOS[i], EPS, 25L),
      numeric(3))
    if (i == 1L) sc <- c(sc, list(V["hill", ]))   # grid control, rho-free
    sc <- c(sc, list(V["br", ]))
  }
  stats <- lapply(sc, function(s) {
    o <- order(s, seq_len(p), na.last = TRUE)
    nul <- s[9:p]
    c(rmax = max(match(1:4, o)),
      snr = (median(nul) - max(s[1:4])) / sd(nul),
      level = median(nul), sep = median(nul) - max(s[1:4]))
  })
  list(model = z$model, n = n, p = p, rho = z$rho,
       stats = do.call(rbind, stats))
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
saveRDS(res, "results/bias_reduced/raw.rds", compress = "xz")

fm <- function(v, k = 3) formatC(v, format = "f", digits = k)
cat("\nHill local a biais reduit de second ordre\n\n")
for (f in cells) {
  out <- res[[f]]; z1 <- out[[1]]
  thr <- qnorm(1 - 1 / (z1$p - 8))
  cat(sprintf("=== %s  n=%d p=%d rho=%.2f, %d replications (seuil de competition %.2f sigma)\n",
              z1$model, z1$n, z1$p, z1$rho, length(out), thr))
  cat(sprintf("%-22s %8s %8s %8s %9s %9s %9s\n", "regle", "Sure-4", "Sure-10",
              "Sure-20", "SNR", "niveau", "ecart"))
  for (i in seq_along(RULES)) {
    g <- function(w) vapply(out, function(z) z$stats[i, w], numeric(1))
    rm_ <- g("rmax")
    cat(sprintf("%-22s %8s %8s %8s %9s %9s %9s\n", RULES[i],
                fm(mean(rm_ <= 4)), fm(mean(rm_ <= 10)), fm(mean(rm_ <= 20)),
                fm(mean(g("snr")), 2), fm(mean(g("level")), 3),
                fm(mean(g("sep")), 4)))
  }
  cat("\n")
}
