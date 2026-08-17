## Yoshida--Umezu conditional-Pickands screening over a tuning grid.
##
## Same statistic as code/R/yoshida_umezu.R, but evaluated over a grid of
## intermediate sequences k at a fixed bandwidth h, with two exact
## optimisations:
##
##   1. the response is sorted ONCE per dataset instead of once per
##      coordinate (the sort does not depend on the coordinate);
##   2. for a given (coordinate, evaluation point, h), the weighted
##      cumulative sum is formed once and all values of k read their three
##      quantiles off it, since only the thresholds k/n, 2k/n, 4k/n depend
##      on k.
##
## Both are algebraic identities, not approximations: the output is
## bitwise identical to looping yu_score_matrix() over the grid.
## tests/test_yu_grid.R verifies this.

source_once <- function(f) if (!exists("pickands_from_quantiles")) source(f)
source_once("code/R/yoshida_umezu.R")

## First index at which a nondecreasing cumulative sum reaches `thr`,
## i.e. which(cw >= thr)[1], vectorised over thr.  `left.open = TRUE`
## makes findInterval count the entries strictly below thr.
first_at_least <- function(cw, thr) findInterval(thr, cw, left.open = TRUE) + 1L

## Scores for one coordinate, all k in `ks`, at bandwidth h.
## `ud` is the coordinate already permuted into decreasing-y order.
yu_scores_coordinate_grid <- function(ud, yd, ks, h, eval_z, marginal) {
  n <- length(yd); nk <- length(ks); nz <- length(eval_z)
  local <- matrix(NA_real_, nz, nk)
  for (iz in seq_len(nz)) {
    zz <- (ud - eval_z[iz]) / h
    w <- (1 - zz^2) * (abs(zz) <= 1)   # Epanechnikov; normalisation cancels
    cw <- cumsum(w); tot <- cw[n]
    if (!is.finite(tot) || tot <= 0) next
    thr <- outer(c(1, 2, 4), ks) / n * tot        # 3 x nk
    at <- first_at_least(cw, as.vector(thr))
    at[at > n] <- NA_integer_
    q <- matrix(yd[at], 3L, nk)
    local[iz, ] <- vapply(seq_len(nk), function(ik)
      pickands_from_quantiles(q[1L, ik], q[2L, ik], q[3L, ik]), 0.0)
  }
  fin <- is.finite(local)
  score <- vapply(seq_len(nk), function(ik) {
    dev <- (local[, ik] - marginal[ik])^2
    mean(dev)                                    # NA if any local is NA
  }, 0.0)
  list(score = score, undefined = 1 - colMeans(fin))
}

## Scores for every coordinate of `u`, all k in `ks`, at bandwidth h.
## Returns a list with p x length(ks) matrices.
yu_score_matrix_grid <- function(u, y, ks, h = 1,
                                 eval_z = seq(.02, .98, length.out = 25)) {
  n <- length(y); p <- ncol(u); ks <- as.integer(ks)
  stopifnot(nrow(u) == n, all(ks >= 1L), all(4L * ks < n), h > 0)
  oy <- order(y, decreasing = TRUE); yd <- y[oy]
  marginal <- vapply(ks, function(k)
    pickands_from_quantiles(yd[k], yd[2L * k], yd[4L * k]), 0.0)
  sc <- matrix(NA_real_, p, length(ks)); un <- sc
  for (j in seq_len(p)) {
    r <- yu_scores_coordinate_grid(u[oy, j], yd, ks, h, eval_z, marginal)
    sc[j, ] <- r$score; un[j, ] <- r$undefined
  }
  dimnames(sc) <- dimnames(un) <- list(NULL, paste0("k", ks))
  list(scores = sc, undefined = un, marginal = marginal)
}
