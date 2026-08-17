## Final comparison campaign, revision 8.
##
## Two changes with respect to campaign6/campaign7.
##
## 1. Selected tuning of the proposed screen moves from (a*, b*) =
##    (0.30, 0.15) to (0.35, 0.15).  The score matrix over the 3x3 block
##    is unchanged, so this is a change of which column is reported, not a
##    change of estimator.
##
## 2. The Yoshida--Umezu competitor is no longer run only at the single
##    tuning transplanted from its own paper.  It is given a grid search
##    of its own (code/R/yu_tuning8.R, run on the SAME datasets used to
##    tune the proposed screen), and the best cell is kept PER MODEL,
##    whereas the proposed screen keeps one pair common to all four
##    models.  Both versions are reported:
##       "Yoshida-Umezu, tuned per model"  best (h, k) for that model
##       "Yoshida-Umezu, paper tuning"     h = 1, k = floor(0.072 n)
##    so that the effect of the grid search on the competitor is visible
##    rather than hidden.
##
## Ranking rules are campaign7's permutation-invariant ones throughout.
##
## Seed stream identical to campaign5/6/7 (131000021 + cell_id*100003 +
## model_index*10007 + r*307), so every replication runs on exactly the
## dataset used by the previous campaigns and the comparison with them is
## paired.  Tie-break stream is campaign7's (973000019 + ...).
##
## usage: Rscript code/R/campaign8.R OUTDIR [MODELS] [NREP] [CORES]
## env:   C8_BLOCK checkpoint block size (default 50)
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "results/campaign8"
MODELS <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else
  c("A1", "A2", "A3", "B1")
NREP  <- if (length(args) >= 3) as.integer(args[3]) else 1000L
CORES <- if (length(args) >= 4) as.integer(args[4]) else 7L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate5.R"); source("code/R/qa_sis.R")
source("code/R/yoshida_umezu.R"); source("code/R/rank_rules.R")

N <- 2000L; RHO <- 0.25; PS <- c(500L, 1000L, 2000L)
ASTAR <- 0.35; BSTAR <- 0.15                       # <- was (0.30, 0.15)
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.10, 0.15, 0.20)
TAUS <- c(0.90, 0.95, 0.99)
MODEL_INDEX <- c(A1 = 1L, B1 = 2L, A2 = 3L, A3 = 4L)
settings <- expand.grid(a = AGRID, b = BGRID)
setl <- data.frame(h = N^(-settings$b) / 2, alpha = N^(-settings$a))
IDX_PUB <- which(settings$a == ASTAR & settings$b == BSTAR)
stopifnot(length(IDX_PUB) == 1L)

## Per-model Yoshida--Umezu tuning, from the grid study.
YU_SEL_FILE <- Sys.getenv("YU8_SELECTED", "results/yu_tuning8/yu_selected.csv")
if (!file.exists(YU_SEL_FILE))
  stop("missing ", YU_SEL_FILE, ": run code/R/yu_tuning8.R first")
yu_sel <- read.csv(YU_SEL_FILE, stringsAsFactors = FALSE)
stopifnot(all(MODELS %in% yu_sel$model))
YU_H <- setNames(yu_sel$h, yu_sel$model)
YU_K <- setNames(as.integer(yu_sel$k), yu_sel$model)
YU_PAPER_H <- 1; YU_PAPER_K <- floor(0.072 * N)
cat("YU per-model tuning in force:\n"); print(yu_sel[, c("model","h","k")])

one_rep <- function(model, p, r, seed) {
  d <- simulate_dataset5(N, p, RHO, model, seed)
  S <- matrix(NA_real_, p, nrow(setl))
  for (k in seq_len(nrow(setl)))
    S[, k] <- vapply(seq_len(p), function(j)
      score_coordinate_cpp(d$z[, j], d$y, setl$h[k], setl$alpha[k],
                           .05)[["score"]], numeric(1))
  cell_id <- match(p, c(2000L, 1000L, 500L))
  tb_seed <- 973000019L + cell_id * 100003L + MODEL_INDEX[[model]] * 10007L +
    r * 307L
  u_tb <- tiebreak_u(p, tb_seed)

  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (N + 1))
  yu_tuned <- yu_score_matrix(uh, d$y, k = YU_K[[model]], h = YU_H[[model]])
  yu_paper <- yu_score_matrix(uh, d$y, k = YU_PAPER_K, h = YU_PAPER_H)

  ords <- c(list(order_selected_new(S, IDX_PUB, u_tb),
                 order_agg_new(S, u_tb),
                 order(-yu_tuned$scores, u_tb, na.last = TRUE),
                 order(-yu_paper$scores, u_tb, na.last = TRUE)),
            lapply(TAUS, function(tt)
              order(-qa_sis_scores(uh, d$y, tau = tt), u_tb)))
  list(model = model, p = p, replicate = r, seed = seed,
       rmax = vapply(ords, function(o) max(match(A_GAMMA, o)), numeric(1)),
       t4g = vapply(ords, function(o) sum(o[1:4] %in% A_GAMMA), numeric(1)),
       t4s = vapply(ords, function(o) sum(o[1:4] %in% A_SCALE), numeric(1)),
       t24s = vapply(ords, function(o) sum(o[1:24] %in% A_SCALE), numeric(1)))
}

RULES <- c("tail selected", "tail aggregated",
           "Yoshida-Umezu tuned", "Yoshida-Umezu paper",
           "quantile .90", "quantile .95", "quantile .99")
BLOCK <- as.integer(Sys.getenv("C8_BLOCK", "50"))
cells <- expand.grid(model = MODELS, p = PS, stringsAsFactors = FALSE)
cells <- cells[order(cells$p), ]
for (ci in seq_len(nrow(cells))) {
  m <- cells$model[ci]; p <- cells$p[ci]
  cell_id <- match(p, c(2000L, 1000L, 500L))
  for (b0 in seq(1L, NREP, by = BLOCK)) {
    b1 <- min(NREP, b0 + BLOCK - 1L)
    path <- file.path(outdir, sprintf("c8_%s_p%04d_r%04d_%04d.rds",
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
files <- list.files(outdir, pattern = "^c8_", full.names = TRUE)
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
cat("CAMPAIGN8 DONE:", length(jobs), "replications summarized\n")
cat("selected pair (a*,b*) = (", ASTAR, ",", BSTAR, ")\n")
