## Draft-5 final comparison campaign on the A1/A2/A3/B1 suite.
##   n = 2000, p in {500, 1000, 2000}, rho = 0.25, 1000 replications.
## Rules (all on the same dataset within a replication):
##   1 tail-index SIS at the selected pair (ASTAR, BSTAR)
##   2 tail-index SIS, minimum-rank aggregation over the 3x3 block N9
##   3 Yoshida-Umezu, transplanted tuning h = 1, k = floor(0.072 n)
##   4-6 quantile SIS at tau = 0.90, 0.95, 0.99
## Metrics: Sure-d (d in 4,10,20,30,50), E(Rmax), Med(Rmax), and the
## top-4 / top-24 composition against A_GAMMA and A_SCALE.
##
## Seed stream identical to campaign5 (131000021 + cell*100003 +
## model_index*10007 + r*307) with the index map A1->1, B1->2, A2->3, A3->4,
## so that the A1 and A2 datasets coincide exactly with campaign5's M1 and
## M3 datasets (verified by tests/test_campaign_equivalence.R).
## usage: Rscript code/R/campaign6.R OUTDIR [MODELS] [NREP] [CORES]
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "results/campaign6"
MODELS <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else
  c("A1", "A2", "A3", "B1")
NREP  <- if (length(args) >= 3) as.integer(args[3]) else 1000L
CORES <- if (length(args) >= 4) as.integer(args[4]) else 7L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate5.R"); source("code/R/qa_sis.R")
source("code/R/yoshida_umezu.R")

N <- 2000L; RHO <- 0.25; PS <- c(500L, 1000L, 2000L)
ASTAR <- 0.30; BSTAR <- 0.15
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.10, 0.15, 0.20)
TAUS <- c(0.90, 0.95, 0.99)
MODEL_INDEX <- c(A1 = 1L, B1 = 2L, A2 = 3L, A3 = 4L)
settings <- expand.grid(a = AGRID, b = BGRID)
setl <- data.frame(h = N^(-settings$b) / 2, alpha = N^(-settings$a))
IDX_PUB <- which(settings$a == ASTAR & settings$b == BSTAR)

one_rep <- function(model, p, r, seed) {
  d <- simulate_dataset5(N, p, RHO, model, seed)
  S <- matrix(NA_real_, p, nrow(setl)); R <- S
  for (k in seq_len(nrow(setl))) {
    S[, k] <- vapply(seq_len(p), function(j)
      score_coordinate_cpp(d$z[, j], d$y, setl$h[k], setl$alpha[k],
                           .05)[["score"]], numeric(1))
    R[, k] <- rank(S[, k], ties.method = "first")
  }
  amin <- apply(R, 1L, min)
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (N + 1))
  yu <- yu_score_matrix(uh, d$y, k = floor(0.072 * N), h = 1)
  ords <- c(list(order(S[, IDX_PUB], seq_len(p), na.last = TRUE),
                 order(amin, seq_len(p)),
                 order(-yu$scores, seq_len(p))),
            lapply(TAUS, function(tt)
              order(-qa_sis_scores(uh, d$y, tau = tt), seq_len(p))))
  list(model = model, p = p, replicate = r, seed = seed,
       rmax = vapply(ords, function(o) max(match(A_GAMMA, o)), numeric(1)),
       t4g = vapply(ords, function(o) sum(o[1:4] %in% A_GAMMA), numeric(1)),
       t4s = vapply(ords, function(o) sum(o[1:4] %in% A_SCALE), numeric(1)),
       t24s = vapply(ords, function(o) sum(o[1:24] %in% A_SCALE), numeric(1)))
}

RULES <- c("tail selected", "tail aggregated", "Yoshida-Umezu",
           "quantile .90", "quantile .95", "quantile .99")
BLOCK <- as.integer(Sys.getenv("C6_BLOCK", "50"))  # set C6_BLOCK >= cores on large machines
cells <- expand.grid(model = MODELS, p = PS, stringsAsFactors = FALSE)
cells <- cells[order(cells$p), ]
for (ci in seq_len(nrow(cells))) {
  m <- cells$model[ci]; p <- cells$p[ci]
  cell_id <- match(p, c(2000L, 1000L, 500L))  # campaign5 cell order (desc. p)
  for (b0 in seq(1L, NREP, by = BLOCK)) {
    b1 <- min(NREP, b0 + BLOCK - 1L)
    path <- file.path(outdir, sprintf("c6_%s_p%04d_r%04d_%04d.rds",
                                      m, p, b0, b1))
    if (file.exists(path)) next
    out <- parallel::mclapply(b0:b1, function(r) {
      seed <- 131000021L + cell_id * 100003L + MODEL_INDEX[[m]] * 10007L +
        r * 307L
      one_rep(m, p, r, seed)
    }, mc.cores = CORES, mc.preschedule = FALSE)
    bad <- vapply(out, function(z) inherits(z, "try-error") || is.null(z),
                  logical(1))
    if (any(bad)) stop("cell ", path, ": ", sum(bad), " failed")
    saveRDS(out, path, compress = "xz")
    cat(format(Sys.time(), "%H:%M:%S"), m, p,
        sprintf("reps %d-%d done\n", b0, b1)); flush.console()
  }
}

## Summary.
files <- list.files(outdir, pattern = "^c6_", full.names = TRUE)
jobs <- unlist(lapply(files, readRDS), recursive = FALSE)
rows <- list()
for (m in unique(vapply(jobs, `[[`, "", "model")))
  for (p in PS) {
    w <- Filter(function(z) z$model == m && z$p == p, jobs)
    if (!length(w)) next
    for (k in seq_along(RULES)) {
      v <- vapply(w, function(z) z$rmax[k], numeric(1))
      rows[[paste(m, p, k)]] <- data.frame(
        n = N, p = p, rho = RHO, model = m, rule = RULES[k], reps = length(w),
        sure4 = mean(v <= 4), sure10 = mean(v <= 10), sure20 = mean(v <= 20),
        sure30 = mean(v <= 30), sure50 = mean(v <= 50),
        ermax = round(mean(v), 3), medmax = median(v),
        top4_gamma = mean(vapply(w, function(z) z$t4g[k], numeric(1))),
        top4_scale = mean(vapply(w, function(z) z$t4s[k], numeric(1))),
        top24_scale = mean(vapply(w, function(z) z$t24s[k], numeric(1))))
    }
  }
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, file.path(outdir, "summary.csv"), row.names = FALSE)
cat("CAMPAIGN6 DONE:", length(jobs), "replications summarized\n")
