## Two-stage pipeline: screen with the quantile screen, then rank the survivors
## with the tail-index screen.
##
## Rationale.  The tail-index score of a coordinate is marginal, so restricting
## attention to a subset does not change any score -- it only changes the
## competition.  The diagnostic in code/R/diag_tuning_snr.R showed that the
## screen fails because an active coordinate must beat the minimum of p-4 null
## scores, which for p = 2000 sits about 3.29 sigma below the null median while
## the realised separation is 3.1-3.5 sigma.  Cutting the field to d1 = 50
## lowers that bar to about 2.05 sigma.  Stage 1 supplies recall (its Sure-50
## is high in every model), stage 2 supplies the tail-index interpretation and
## demotes the scale-only variables stage 1 lets through.
##
## The pipeline can only lose what stage 1 discards, so P(all actives survive
## stage 1) is reported as the ceiling on its Sure-d.
##
## Replications are replayed from the seeds stored in the grid cells, so the
## baselines here must reproduce the campaign's ranks exactly; that is asserted.
##
## usage: KAPPA=0.20 Rscript code/R/test_pipeline.R CORES CELL.rds [CELL.rds ...]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R")
source("code/R/qa_sis.R")

args <- commandArgs(trailingOnly = TRUE)
CORES <- as.integer(args[1L])
cells <- args[-1L]
EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
TAU <- 0.95
D1 <- c(25L, 50L, 100L, 200L)      # stage-1 set sizes
DS <- c(4, 10, 20)

one <- function(z) {
  n <- z$n; p <- z$p
  d <- simulate_dataset3(n, p, z$rho, z$model, z$seed)
  h <- n^(-BSTAR) / 2; alpha <- n^(-ASTAR)
  sct <- vapply(seq_len(p), function(j)
    score_coordinate_cpp(d$z[, j], d$y, h, alpha, EPS), numeric(5))["score", ]
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (n + 1))
  scq <- qa_sis_scores(uh, d$y, tau = TAU)
  ord_t <- order(sct, seq_len(p), na.last = TRUE)
  ord_q <- order(-scq, seq_len(p))
  ## the replay must reproduce the campaign exactly
  stopifnot(identical(as.integer(match(1:4, ord_t)), as.integer(z$ranks_ours)),
            identical(as.integer(match(1:4, ord_q)),
                      as.integer(z$q$q950$ranks)))
  pipe <- vapply(D1, function(d1) {
    surv <- ord_q[seq_len(d1)]
    if (!all(1:4 %in% surv)) return(Inf)
    ord2 <- surv[order(sct[surv], surv)]
    max(match(1:4, ord2))
  }, numeric(1))
  survived <- vapply(D1, function(d1) all(1:4 %in% ord_q[seq_len(d1)]),
                     logical(1))
  list(model = z$model, n = n, p = p, rho = z$rho,
       rmax_t = max(match(1:4, ord_t)), rmax_q = max(match(1:4, ord_q)),
       rmax_pipe = pipe, survived = survived)
}

res <- list()
for (f in cells) {
  x <- readRDS(f)
  out <- parallel::mclapply(x, one, mc.cores = CORES, mc.preschedule = FALSE)
  bad <- vapply(out, function(z) inherits(z, "try-error") || is.null(z),
                logical(1))
  if (any(bad)) stop(f, ": ", sum(bad), " replicates failed: ",
                     as.character(out[bad][[1]]))
  res[[f]] <- out
  cat(format(Sys.time(), "%H:%M:%S"), basename(f), "done\n"); flush.console()
}
saveRDS(res, "results/pipeline/raw.rds", compress = "xz")

fm <- function(v) formatC(v, format = "f", digits = 3)
cat(sprintf("\nPipeline : quantile SIS (tau=%.2f) en etage 1, puis tail-index SIS\n\n", TAU))
for (f in cells) {
  out <- res[[f]]; z1 <- out[[1]]
  cat(sprintf("=== %s  n=%d p=%d rho=%.2f, %d replications\n",
              z1$model, z1$n, z1$p, z1$rho, length(out)))
  cat(sprintf("%-22s %8s %8s %8s %10s\n", "regle", "Sure-4", "Sure-10",
              "Sure-20", "survie e1"))
  rt <- vapply(out, `[[`, numeric(1), "rmax_t")
  rq <- vapply(out, `[[`, numeric(1), "rmax_q")
  cat(sprintf("%-22s %8s %8s %8s %10s\n", "tail-index seul",
              fm(mean(rt <= 4)), fm(mean(rt <= 10)), fm(mean(rt <= 20)), "-"))
  cat(sprintf("%-22s %8s %8s %8s %10s\n", "quantile seul",
              fm(mean(rq <= 4)), fm(mean(rq <= 10)), fm(mean(rq <= 20)), "-"))
  for (i in seq_along(D1)) {
    rp <- vapply(out, function(z) z$rmax_pipe[i], numeric(1))
    sv <- mean(vapply(out, function(z) z$survived[i], logical(1)))
    cat(sprintf("%-22s %8s %8s %8s %10s\n",
                sprintf("pipeline d1=%d", D1[i]),
                fm(mean(rp <= 4)), fm(mean(rp <= 10)), fm(mean(rp <= 20)),
                fm(sv)))
  }
  cat("\n")
}
