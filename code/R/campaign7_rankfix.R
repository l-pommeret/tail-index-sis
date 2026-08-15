## Rank-fix campaign: paired before/after comparison of the tie-break rule
## for the tail-index screens, on the exact datasets of campaign5/campaign6
## (same seed stream, verified per-replication against the stored results).
##
## Per replication, from ONE 9-setting score matrix:
##   old_sel / old_agg  exact campaign6 rules (column-index tie-break)
##   new_sel / new_agg  permutation-invariant rules (lexicographic key
##                      min/median/mean/max of the 9 ranks + seeded random
##                      residual tie-break; see code/R/rank_rules.R)
##   prm_sel / prm_agg  new rules after a seeded random permutation of the
##                      columns, metrics evaluated against the permuted
##                      positions of the actives.  Because every screen is
##                      marginal, permuting the rows of the score matrix is
##                      bitwise identical to rerunning on permuted columns
##                      (verified by tests/test_rank_invariance.R); this
##                      arm therefore costs no extra scoring pass.
## The four competitor rules (Yoshida-Umezu, quantile .90/.95/.99) are NOT
## recomputed: their scores are continuous (no exact ties: see
## results/diag_ties6/competitor_ties.csv) and marginal, so both the
## tie-break change and the column permutation leave their results
## unchanged; the stored campaign5/campaign6 numbers remain valid.
##
## Dedicated seed streams (never used elsewhere):
##   tie-break    973000019 + cell_id*100003 + model_index*10007 + r*307
##   permutation  987000037 + cell_id*100003 + model_index*10007 + r*307
## Dataset seeds are campaign5's: 131000021 + cell_id*100003 +
## model_index*10007 + r*307, cell_id = match(p, c(2000,1000,500)),
## model index A1->1, B1->2, A2->3, A3->4.
##
## usage: Rscript code/R/campaign7_rankfix.R OUTDIR [MODELS] [NREP] [CORES]
## env:   C7_BLOCK checkpoint block size (default 50; >= cores on the
##        cluster, e.g. C7_BLOCK=250)
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "results/campaign7"
MODELS <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else
  c("A1", "A2", "A3", "B1")
NREP  <- if (length(args) >= 3) as.integer(args[3]) else 1000L
CORES <- if (length(args) >= 4) as.integer(args[4]) else 7L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate5.R"); source("code/R/rank_rules.R")

N <- 2000L; RHO <- 0.25; PS <- c(500L, 1000L, 2000L)
ASTAR <- 0.30; BSTAR <- 0.15
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.10, 0.15, 0.20)
MODEL_INDEX <- c(A1 = 1L, B1 = 2L, A2 = 3L, A3 = 4L)
settings <- expand.grid(a = AGRID, b = BGRID)
setl <- data.frame(h = N^(-settings$b) / 2, alpha = N^(-settings$a))
IDX_PUB <- which(settings$a == ASTAR & settings$b == BSTAR)
SEED_BASE <- 131000021L; TB_BASE <- 973000019L; PM_BASE <- 987000037L

perm_of <- function(p, seed) {
  old <- if (exists(".Random.seed", .GlobalEnv))
    get(".Random.seed", .GlobalEnv) else NULL
  set.seed(seed); s <- sample.int(p)
  if (is.null(old)) rm(".Random.seed", envir = .GlobalEnv) else
    assign(".Random.seed", old, .GlobalEnv)
  s
}

one_rep <- function(model, p, r) {
  cell_id <- match(p, c(2000L, 1000L, 500L))
  off <- cell_id * 100003L + MODEL_INDEX[[model]] * 10007L + r * 307L
  d <- simulate_dataset5(N, p, RHO, model, SEED_BASE + off)
  S <- matrix(NA_real_, p, nrow(setl))
  for (k in seq_len(nrow(setl)))
    S[, k] <- vapply(seq_len(p), function(j)
      score_coordinate_cpp(d$z[, j], d$y, setl$h[k], setl$alpha[k],
                           .05)[["score"]], numeric(1))
  u <- tiebreak_u(p, TB_BASE + off)
  ords <- list(old_sel = order_selected_old(S, IDX_PUB),
               old_agg = order_agg_old(S),
               new_sel = order_selected_new(S, IDX_PUB, u),
               new_agg = order_agg_new(S, u))
  M <- vapply(ords, metrics_from_order, numeric(4),
              actives = A_GAMMA, scale_set = A_SCALE)
  perm <- perm_of(p, PM_BASE + off)
  Sp <- S[perm, , drop = FALSE]
  Ag <- match(A_GAMMA, perm); As <- match(A_SCALE, perm)
  Mp <- vapply(list(prm_sel = order_selected_new(Sp, IDX_PUB, u),
                    prm_agg = order_agg_new(Sp, u)),
               metrics_from_order, numeric(4), actives = Ag, scale_set = As)
  list(model = model, p = p, replicate = r,
       metrics = cbind(M, Mp), ties = tie_stats(S, A_GAMMA))
}

BLOCK <- as.integer(Sys.getenv("C7_BLOCK", "50"))
cells <- expand.grid(model = MODELS, p = PS, stringsAsFactors = FALSE)
cells <- cells[order(cells$p), ]
for (ci in seq_len(nrow(cells))) {
  m <- cells$model[ci]; p <- cells$p[ci]
  for (b0 in seq(1L, NREP, by = BLOCK)) {
    b1 <- min(NREP, b0 + BLOCK - 1L)
    path <- file.path(outdir, sprintf("c7_%s_p%04d_r%04d_%04d.rds",
                                      m, p, b0, b1))
    if (file.exists(path)) next
    out <- parallel::mclapply(b0:b1, function(r) one_rep(m, p, r),
                              mc.cores = CORES, mc.preschedule = FALSE)
    bad <- vapply(out, function(z) inherits(z, "try-error") || is.null(z),
                  logical(1))
    if (any(bad)) stop("cell ", path, ": ", sum(bad), " failed")
    saveRDS(out, path, compress = "xz")
    cat(format(Sys.time(), "%H:%M:%S"), m, p,
        sprintf("reps %d-%d done\n", b0, b1)); flush.console()
  }
}

## ---- summary + paired before/after table -------------------------------
files <- list.files(outdir, pattern = "^c7_", full.names = TRUE)
jobs <- unlist(lapply(files, readRDS), recursive = FALSE)
ARMS <- c("old_sel", "old_agg", "new_sel", "new_agg", "prm_sel", "prm_agg")
srows <- list(); brows <- list(); crows <- list()
for (m in unique(vapply(jobs, `[[`, "", "model"))) for (p in PS) {
  w <- Filter(function(z) z$model == m && z$p == p, jobs)
  if (!length(w)) next
  rmx <- sapply(w, function(z) z$metrics["rmax", ])      # ARMS x reps
  for (a in ARMS) {
    v <- rmx[a, ]
    srows[[paste(m, p, a)]] <- data.frame(
      n = N, p = p, model = m, arm = a, reps = length(w),
      sure4 = mean(v <= 4), sure10 = mean(v <= 10), sure20 = mean(v <= 20),
      sure30 = mean(v <= 30), sure50 = mean(v <= 50),
      ermax = round(mean(v), 3), medmax = median(v),
      top4_gamma = mean(sapply(w, function(z) z$metrics["t4g", a])),
      top4_scale = mean(sapply(w, function(z) z$metrics["t4s", a])),
      top24_scale = mean(sapply(w, function(z) z$metrics["t24s", a])))
  }
  pse <- function(x, y) sd(x - y) / sqrt(length(x))     # paired SE
  for (cmp in list(c("old_agg", "new_agg"), c("old_agg", "prm_agg"),
                   c("old_sel", "new_sel"))) {
    vo <- rmx[cmp[1], ]; vn <- rmx[cmp[2], ]
    brows[[paste(m, p, cmp[2])]] <- data.frame(
      model = m, p = p, comparison = paste(cmp[1], "vs", cmp[2]),
      reps = length(vo),
      sure4_old = mean(vo <= 4), sure4_new = mean(vn <= 4),
      d_sure4 = mean((vn <= 4) - (vo <= 4)),
      se_sure4 = pse(vn <= 4, vo <= 4),
      sure20_old = mean(vo <= 20), sure20_new = mean(vn <= 20),
      d_sure20 = mean((vn <= 20) - (vo <= 20)),
      se_sure20 = pse(vn <= 20, vo <= 20),
      ermax_old = round(mean(vo), 3), ermax_new = round(mean(vn), 3),
      d_ermax = round(mean(vn - vo), 3), se_ermax = round(pse(vn, vo), 3),
      medmax_old = median(vo), medmax_new = median(vn),
      changed_reps = mean(vn != vo))
  }
  ties <- sapply(w, `[[`, "ties")
  crows[[paste(m, p)]] <- data.frame(
    model = m, p = p, t(rowMeans(ties)))
}
write.csv(do.call(rbind, srows), file.path(outdir, "summary.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, brows), file.path(outdir, "before_after.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, crows), file.path(outdir, "tie_stats.csv"),
          row.names = FALSE)

## ---- cross-check: old rules must reproduce the stored campaigns --------
chk <- function(stored, model_lab, model_new, src) {
  for (p in PS) {
    w5 <- Filter(function(z) z$model == model_lab && z$p == p, stored)
    w7 <- Filter(function(z) z$model == model_new && z$p == p, jobs)
    if (!length(w5) || !length(w7)) next
    r5 <- vapply(w5, `[[`, 0L, "replicate")
    m5 <- sapply(w5, function(z) z$rmax[1:2])            # rules 1-2
    r7 <- vapply(w7, `[[`, 0L, "replicate")
    m7 <- sapply(w7, function(z) z$metrics["rmax", c("old_sel", "old_agg")])
    common <- intersect(r5, r7)
    dd <- max(abs(m5[, match(common, r5)] - m7[, match(common, r7)]))
    cat(sprintf("crosscheck %s p=%d vs %s: %d common reps, max|diff| = %g\n",
                model_new, p, src, length(common), dd))
    if (dd != 0) warning("old-rule metrics do not reproduce ", src)
  }
}
try({
  f6 <- list.files("results/campaign6", pattern = "^c6_", full.names = TRUE)
  if (length(f6)) {
    j6 <- unlist(lapply(f6, readRDS), recursive = FALSE)
    for (m in intersect(MODELS, c("A3", "B1"))) chk(j6, m, m, "campaign6")
  }
  f5 <- list.files("results/campaign5", pattern = "^c5_", full.names = TRUE)
  if (length(f5)) {
    j5 <- unlist(lapply(f5, readRDS), recursive = FALSE)
    if ("A1" %in% MODELS) chk(j5, "M1", "A1", "campaign5")
    if ("A2" %in% MODELS) chk(j5, "M3", "A2", "campaign5")
  }
})
cat("CAMPAIGN7 DONE:", length(jobs), "replications summarized\n")
