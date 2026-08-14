## Evaluate the extremal quantile-regression slope screen against the two
## existing screens, on replications replayed from the grid cells' seeds.
## usage: KAPPA=0.20 Rscript code/R/test_qr_slope.R CORES CELL.rds [CELL.rds ...]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R")
source("code/R/qa_sis.R")
source("code/R/qr_slope.R")

args <- commandArgs(trailingOnly = TRUE)
CORES <- as.integer(args[1L])
cells <- args[-1L]
dir.create("results/qr_slope", showWarnings = FALSE, recursive = TRUE)
EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
TAUS <- c(0.90, 0.95, 0.975, 0.99)
DS <- c(4, 10, 20)
RULES <- c("tail-index", "quantile .95", "QR slope (moyenne)", "QR slope (sd)")

one <- function(z) {
  n <- z$n; p <- z$p
  d <- simulate_dataset3(n, p, z$rho, z$model, z$seed)
  h <- n^(-BSTAR) / 2; alpha <- n^(-ASTAR)
  sct <- vapply(seq_len(p), function(j)
    score_coordinate_cpp(d$z[, j], d$y, h, alpha, EPS), numeric(5))["score", ]
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (n + 1))
  scq <- qa_sis_scores(uh, d$y, tau = 0.95)
  qs <- qr_slope_scores(uh, d$y, taus = TAUS)
  ords <- list(
    order(sct, seq_len(p), na.last = TRUE),                  # small = active
    order(-scq, seq_len(p)),                                 # large = active
    order(qs[, "slope_mean"], seq_len(p), na.last = TRUE),   # small = active
    order(-qs[, "slope_sd"], seq_len(p), na.last = TRUE))    # large = active
  list(model = z$model, n = n, p = p, rho = z$rho,
       rmax = vapply(ords, function(o) max(match(1:4, o)), numeric(1)),
       t4g = vapply(ords, function(o) sum(o[1:4] %in% 1:4), numeric(1)),
       t4s = vapply(ords, function(o) sum(o[1:4] %in% 5:8), numeric(1)),
       na_slope = mean(is.na(qs[, "slope_mean"])))
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
dir.create("results/qr_slope", showWarnings = FALSE, recursive = TRUE)
saveRDS(res, "results/qr_slope/raw.rds", compress = "xz")

fm <- function(v) formatC(v, format = "f", digits = 3)
cat("\nScreen par pente de regression quantile extreme\n")
cat(sprintf("tau = %s, B-splines df=3, sur log(Y)\n\n",
            paste(TAUS, collapse = ", ")))
for (f in cells) {
  out <- res[[f]]; z1 <- out[[1]]
  cat(sprintf("=== %s  n=%d p=%d rho=%.2f, %d replications (pentes NA : %.3f)\n",
              z1$model, z1$n, z1$p, z1$rho, length(out),
              mean(vapply(out, `[[`, numeric(1), "na_slope"))))
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
