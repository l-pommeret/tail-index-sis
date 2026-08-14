## Evaluate the exceedance-regression screen against the two existing screens,
## on replications replayed from the grid cells' seeds.
## usage: KAPPA=0.20 Rscript code/R/test_exceedance.R CORES CELL.rds [CELL.rds ...]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R")
source("code/R/qa_sis.R")
source("code/R/exceedance_hill.R")

args <- commandArgs(trailingOnly = TRUE)
CORES <- as.integer(args[1L])
cells <- args[-1L]
dir.create("results/exceedance", showWarnings = FALSE, recursive = TRUE)
EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
TAU0 <- c(0.90, 0.95)
RULES <- c("tail-index", "quantile .95",
           paste0("exc t0=", rep(TAU0, each = 3), c(" (moy)", " (sd)", " (F)")))

one <- function(z) {
  n <- z$n; p <- z$p
  d <- simulate_dataset3(n, p, z$rho, z$model, z$seed)
  h <- n^(-BSTAR) / 2; alpha <- n^(-ASTAR)
  sct <- vapply(seq_len(p), function(j)
    score_coordinate_cpp(d$z[, j], d$y, h, alpha, EPS), numeric(5))["score", ]
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (n + 1))
  scq <- qa_sis_scores(uh, d$y, tau = 0.95)
  ords <- list(order(sct, seq_len(p), na.last = TRUE),
               order(-scq, seq_len(p)))
  for (t0 in TAU0) {
    e <- exceedance_hill_scores(uh, d$y, tau0 = t0)
    ords <- c(ords, list(
      order(e[, "g_mean"], seq_len(p), na.last = TRUE),
      order(-e[, "g_sd"], seq_len(p), na.last = TRUE),
      order(-e[, "g_F"], seq_len(p), na.last = TRUE)))
  }
  list(model = z$model, n = n, p = p, rho = z$rho,
       rmax = vapply(ords, function(o) max(match(1:4, o)), numeric(1)),
       t4g = vapply(ords, function(o) sum(o[1:4] %in% 1:4), numeric(1)),
       t4s = vapply(ords, function(o) sum(o[1:4] %in% 5:8), numeric(1)))
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
saveRDS(res, "results/exceedance/raw.rds", compress = "xz")

fm <- function(v) formatC(v, format = "f", digits = 3)
cat("\nScreen par regression sur les depassements\n\n")
for (f in cells) {
  out <- res[[f]]; z1 <- out[[1]]
  cat(sprintf("=== %s  n=%d p=%d rho=%.2f, %d replications\n",
              z1$model, z1$n, z1$p, z1$rho, length(out)))
  cat(sprintf("%-20s %8s %8s %8s %10s %10s\n", "regle", "Sure-4", "Sure-10",
              "Sure-20", "top4 gamma", "top4 ech."))
  for (i in seq_along(RULES)) {
    rm_ <- vapply(out, function(z) z$rmax[i], numeric(1))
    cat(sprintf("%-20s %8s %8s %8s %10s %10s\n", RULES[i],
                fm(mean(rm_ <= 4)), fm(mean(rm_ <= 10)), fm(mean(rm_ <= 20)),
                fm(mean(vapply(out, function(z) z$t4g[i], numeric(1)))),
                fm(mean(vapply(out, function(z) z$t4s[i], numeric(1))))))
  }
  cat("\n")
}
