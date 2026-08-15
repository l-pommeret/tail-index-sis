## Permutation-invariance test for the ranking rules (audit deliverable).
##
## For a fixed dataset (X, Y):
##   1. run the full scoring + ranking pipeline on the original columns;
##   2. permute the columns of X, rerun the identical pipeline, map the
##      resulting ordering back to the original labels;
##   3. the NEW rule must give the same ranking up to reproducible
##      randomness restricted to exact full-key ties;
##   4. the OLD rule must be shown to depend on the column labels whenever
##      a mixed tie exists (the test FAILS if the index dependence has any
##      effect and the invariant rule does not remove it).
## Also verifies bitwise that every screen's score vector is equivariant
## under column permutation (all screens are marginal).
suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate5.R"); source("code/R/rank_rules.R")

N <- 2000L; RHO <- 0.25; P <- 500L
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.10, 0.15, 0.20)
settings <- expand.grid(a = AGRID, b = BGRID)
setl <- data.frame(h = N^(-settings$b) / 2, alpha = N^(-settings$a))

score_matrix <- function(z, y) {
  S <- matrix(NA_real_, ncol(z), nrow(setl))
  for (k in seq_len(nrow(setl)))
    S[, k] <- vapply(seq_len(ncol(z)), function(j)
      score_coordinate_cpp(z[, j], y, setl$h[k], setl$alpha[k],
                           .05)[["score"]], numeric(1))
  S
}

## Rank vector (position of each original coordinate) from an ordering of
## permuted columns: column i of Xp is original column perm[i].
rank_of_original <- function(o, perm) {
  pos <- integer(length(o)); pos[o] <- seq_along(o)   # position of Xp-col i
  pos_orig <- pos[match(seq_along(o), perm)]          # position of orig j
  pos_orig
}

set.seed(20260815)  # test-level seed: dataset choice and permutations
n_datasets <- 3L; n_perms <- 2L
old_dependence_seen <- FALSE
for (dsi in seq_len(n_datasets)) {
  model <- c("A1", "A2", "B1")[dsi]
  d <- simulate_dataset5(N, P, RHO, model, seed = 555000011L + dsi)
  S <- score_matrix(d$z, d$y)
  key <- agg_key(setting_ranks_avg(S))
  u <- tiebreak_u(P, 973000019L + dsi)
  o_new <- order_agg_new(S, u)
  r_new <- integer(P); r_new[o_new] <- seq_len(P)
  o_old <- order_agg_old(S)
  r_old <- integer(P); r_old[o_old] <- seq_len(P)
  for (pmi in seq_len(n_perms)) {
    perm <- sample.int(P)
    Sp <- score_matrix(d$z[, perm], d$y)
    ## (a) marginal equivariance, bitwise
    stopifnot(identical(Sp, S[perm, ]))
    ## (b) NEW rule: same tie-break seed, ranking mapped back must agree
    ## except within exact full-key tie groups.
    up <- tiebreak_u(P, 973000019L + dsi)
    r_new_p <- rank_of_original(order_agg_new(Sp, up), perm)
    disagree <- which(r_new_p != r_new)
    if (length(disagree)) {
      keystr <- apply(key, 1L, paste, collapse = "|")
      grp_ok <- all(vapply(disagree, function(j)
        sum(keystr == keystr[j]) > 1L, logical(1)))
      stopifnot(grp_ok)
      ## and the multiset of positions within each tie group is preserved
      for (ks in unique(keystr[disagree])) {
        g <- which(keystr == ks)
        stopifnot(identical(sort(r_new[g]), sort(r_new_p[g])))
      }
    }
    ## (c) OLD rule: record whether the mapped-back ranking changes.
    r_old_p <- rank_of_original(order_agg_old(Sp), perm)
    if (!identical(r_old_p, r_old)) old_dependence_seen <- TRUE
  }
  cat(sprintf("PASS %s: new rule invariant over %d permutations\n",
              model, n_perms))
}
if (!old_dependence_seen)
  cat("NOTE: no amin tie was hit by these permutations; old-rule index\n",
      "dependence not exercised on these datasets\n")
cat(sprintf("OLD rule label-dependent on these datasets: %s\n",
            old_dependence_seen))
cat("ALL RANK-INVARIANCE CHECKS PASS\n")
