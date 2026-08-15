## Tie diagnostic for the aggregated min-rank screen (audit step, small
## replication count): quantifies how often the OLD column-index tie-break
## can move the reported metrics, before any correction is adopted.
## Datasets and seeds are exactly those of campaign5/campaign6.
## usage: Rscript code/R/diag_ties6.R OUTDIR [NREP] [CORES] [COMPET_REPS]
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "results/diag_ties6"
NREP   <- if (length(args) >= 2) as.integer(args[2]) else 25L
CORES  <- if (length(args) >= 3) as.integer(args[3]) else 7L
COMPET <- if (length(args) >= 4) as.integer(args[4]) else 5L
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
TIEBREAK_BASE <- 973000019L  # dedicated stream, never used elsewhere

score_matrix <- function(d, p) {
  S <- matrix(NA_real_, p, nrow(setl))
  for (k in seq_len(nrow(setl)))
    S[, k] <- vapply(seq_len(p), function(j)
      score_coordinate_cpp(d$z[, j], d$y, setl$h[k], setl$alpha[k],
                           .05)[["score"]], numeric(1))
  S
}

one <- function(model, p, r) {
  cell_id <- match(p, c(2000L, 1000L, 500L))
  seed <- 131000021L + cell_id * 100003L + MODEL_INDEX[[model]] * 10007L +
    r * 307L
  tb <- TIEBREAK_BASE + cell_id * 100003L + MODEL_INDEX[[model]] * 10007L +
    r * 307L
  d <- simulate_dataset5(N, p, RHO, model, seed)
  S <- score_matrix(d, p)
  u <- tiebreak_u(p, tb)
  ords <- list(old  = order_agg_old(S),
               fav  = order_agg_bound(S, A_GAMMA, TRUE),
               adv  = order_agg_bound(S, A_GAMMA, FALSE),
               rnd  = order_agg_minrandom(S, u),
               new  = order_agg_new(S, u))
  M <- vapply(ords, metrics_from_order, numeric(4),
              actives = A_GAMMA, scale_set = A_SCALE)
  sel_old <- metrics_from_order(order_selected_old(S, IDX_PUB),
                                A_GAMMA, A_SCALE)
  sel_new <- metrics_from_order(order_selected_new(S, IDX_PUB, u),
                                A_GAMMA, A_SCALE)
  ts <- tie_stats(S, A_GAMMA)
  data.frame(model = model, p = p, replicate = r, t(ts),
             t(setNames(as.vector(M), paste(rep(colnames(M), each = 4),
                                            rownames(M), sep = "_"))),
             sel_old_rmax = sel_old[["rmax"]], sel_new_rmax = sel_new[["rmax"]],
             sel_ties_differ = !identical(sel_old, sel_new))
}

cells <- expand.grid(model = names(MODEL_INDEX), p = PS,
                     stringsAsFactors = FALSE)
rows <- list()
for (ci in seq_len(nrow(cells))) {
  m <- cells$model[ci]; p <- cells$p[ci]
  out <- parallel::mclapply(seq_len(NREP), function(r) one(m, p, r),
                            mc.cores = CORES, mc.preschedule = FALSE)
  rows[[ci]] <- do.call(rbind, out)
  cat(format(Sys.time(), "%H:%M:%S"), m, p, "done\n"); flush.console()
}
diag <- do.call(rbind, rows)
write.csv(diag, file.path(outdir, "diag.csv"), row.names = FALSE)

## Competitor tie check (their scores are continuous; exact ties should be
## absent, in which case the index tie-break is inert for them).
source("code/R/qa_sis.R"); source("code/R/yoshida_umezu.R")
comp <- do.call(rbind, parallel::mclapply(seq_len(COMPET), function(r) {
  p <- 1000L; cell_id <- 2L
  do.call(rbind, lapply(names(MODEL_INDEX), function(m) {
    seed <- 131000021L + cell_id * 100003L + MODEL_INDEX[[m]] * 10007L +
      r * 307L
    d <- simulate_dataset5(N, p, RHO, m, seed)
    uh <- apply(d$z, 2L, function(x)
      rank(x, ties.method = "average") / (N + 1))
    yu <- yu_score_matrix(uh, d$y, k = floor(0.072 * N), h = 1)$scores
    qs <- sapply(c(.90, .95, .99), function(tt)
      qa_sis_scores(uh, d$y, tau = tt))
    data.frame(model = m, replicate = r,
               yu_dup = sum(duplicated(yu)), yu_na = sum(!is.finite(yu)),
               q_dup = sum(apply(qs, 2, function(v) sum(duplicated(v)))),
               q_na = sum(!is.finite(qs)))
  }))
}, mc.cores = CORES))
write.csv(comp, file.path(outdir, "competitor_ties.csv"), row.names = FALSE)

## Printed summary.
agg <- aggregate(cbind(exact_score_ties, n_inf, tied_coords_amin,
                       mixed4, mixed20,
                       top4_change_adv = old_t4g != adv_t4g,
                       top4_change_new = old_t4g != new_t4g,
                       top20_change_adv = (old_rmax <= 20) != (adv_rmax <= 20),
                       top20_change_new = (old_rmax <= 20) != (new_rmax <= 20),
                       rmax_change_adv = old_rmax != adv_rmax,
                       rmax_change_new = old_rmax != new_rmax,
                       old_is_favourable = old_rmax == fav_rmax &
                         old_t4g == fav_t4g) ~ model + p, diag, mean)
print(agg, digits = 3)
cat("\nSelected-pair rule: orderings differ old vs new in",
    sum(diag$sel_ties_differ), "of", nrow(diag), "replications\n")
cat("Competitor exact ties (should be 0):\n")
print(colSums(comp[, 3:6]))
cat("WROTE", file.path(outdir, "diag.csv"), "\n")
