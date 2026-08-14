## Best of the two working ideas combined: stage 1 screens with the quantile
## screen, stage 2 ranks the survivors by the MEDIAN RANK of the tail-index
## score over a 3x3 tuning grid rather than by a single tuning.  Both ingredients
## improved things on their own; this asks whether they compose.
## usage: KAPPA=0.20 Rscript code/R/test_pipeline_agg.R CORES CELL.rds [...]
suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R"); source("code/R/qa_sis.R")
args <- commandArgs(trailingOnly = TRUE)
CORES <- as.integer(args[1L]); cells <- args[-1L]
dir.create("results/pipeline_agg", showWarnings = FALSE, recursive = TRUE)
EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.05, 0.10, 0.15)
D1 <- c(25L, 50L)
score_vec <- function(z, y, a, b) {
  n <- nrow(z); p <- ncol(z); h <- n^(-b)/2; alpha <- n^(-a)
  vapply(seq_len(p), function(j)
    score_coordinate_cpp(z[, j], y, h, alpha, EPS), numeric(5))["score", ]
}
one <- function(z0) {
  n <- z0$n; p <- z0$p
  d <- simulate_dataset3(n, p, z0$rho, z0$model, z0$seed)
  base <- score_vec(d$z, d$y, ASTAR, BSTAR)
  R <- matrix(NA_real_, p, length(AGRID) * length(BGRID)); ii <- 0L
  for (a in AGRID) for (b in BGRID) {
    ii <- ii + 1L
    R[, ii] <- rank(score_vec(d$z, d$y, a, b), ties.method = "first",
                    na.last = TRUE)
  }
  agg <- apply(R, 1L, median)
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average")/(n + 1))
  ord_q <- order(-qa_sis_scores(uh, d$y, tau = 0.95), seq_len(p))
  out <- c(max(match(1:4, order(base, seq_len(p), na.last = TRUE))),
           max(match(1:4, order(agg, seq_len(p)))),
           max(match(1:4, ord_q)))
  for (d1 in D1) {
    surv <- ord_q[seq_len(d1)]
    for (sc in list(base, agg)) {
      out <- c(out, if (!all(1:4 %in% surv)) Inf
                    else max(match(1:4, surv[order(sc[surv], surv)])))
    }
  }
  list(model = z0$model, n = n, p = p, rho = z0$rho, rmax = out)
}
RULES <- c("tail-index seul", "tail-index 9 reglages", "quantile .95",
           unlist(lapply(D1, function(d) paste0("pipeline d1=", d,
                                                c("", " + 9 reglages")))))
res <- list()
for (f in cells) {
  x <- readRDS(f)
  o <- parallel::mclapply(x, one, mc.cores = CORES, mc.preschedule = FALSE)
  bad <- vapply(o, function(z) inherits(z, "try-error") || is.null(z), logical(1))
  if (any(bad)) stop(f, ": ", sum(bad), " failed: ", as.character(o[bad][[1]]))
  res[[f]] <- o
  cat(format(Sys.time(), "%H:%M:%S"), basename(f), "done\n"); flush.console()
}
saveRDS(res, "results/pipeline_agg/raw.rds", compress = "xz")
fm <- function(v) formatC(v, format = "f", digits = 3)
cat("\nPipeline et agregation de reglages combines\n\n")
for (f in cells) {
  o <- res[[f]]; z1 <- o[[1]]
  cat(sprintf("=== %s  n=%d p=%d rho=%.2f, %d replications\n", z1$model, z1$n,
              z1$p, z1$rho, length(o)))
  cat(sprintf("%-28s %8s %8s %8s\n", "regle", "Sure-4", "Sure-10", "Sure-20"))
  for (i in seq_along(RULES)) {
    v <- vapply(o, function(z) z$rmax[i], numeric(1))
    cat(sprintf("%-28s %8s %8s %8s\n", RULES[i], fm(mean(v <= 4)),
                fm(mean(v <= 10)), fm(mean(v <= 20))))
  }
  cat("\n")
}
