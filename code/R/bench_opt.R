## Candidate optimizations for the comparison drivers, checked for exact
## equality against the reference implementations before any timing claim.
##
##  (1) Yoshida--Umezu: hoist the response sort out of the coordinate loop
##      (it does not depend on j) and compute the weight cumsum once per
##      evaluation point instead of once per tail probability.
##  (2) Quantile SIS: build the cubic B-spline basis once per coordinate and
##      reuse it across the four values of tau, instead of rebuilding it for
##      each tau.  Same basis, same rq.fit.br calls.
##
## usage: Rscript code/R/bench_opt.R [n] [p]

suppressPackageStartupMessages({ library(quantreg); library(splines) })
source("code/R/generate.R")
source("code/R/yoshida_umezu.R")
source("code/R/qa_sis.R")

args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) >= 1L) as.integer(args[1L]) else 2000L
p <- if (length(args) >= 2L) as.integer(args[2L]) else 200L
TAUS <- c(0.90, 0.95, 0.975, 0.99)

## ---------------------------------------------------------------- (1) YU ----
yu_score_matrix_fast <- function(u, y, k, h = 1,
                                 eval_z = seq(.02, .98, length.out = 25)) {
  n <- length(y)
  stopifnot(nrow(u) == n, k >= 1L, 4L * k < n, h > 0)
  oy <- order(y, decreasing = TRUE)      # identical for every coordinate
  yd <- y[oy]
  probs <- c(k, 2 * k, 4 * k) / n
  marginal <- pickands_from_quantiles(yd[k], yd[2 * k], yd[4 * k])
  out <- vapply(seq_len(ncol(u)), function(j) {
    ud <- u[oy, j]
    local <- vapply(eval_z, function(z) {
      zz <- (ud - z) / h
      w <- (1 - zz^2) * (abs(zz) <= 1)   # Epanechnikov, normalization cancels
      cw <- cumsum(w); total <- cw[n]    # one pass instead of three
      if (!is.finite(total) || total <= 0) return(NA_real_)
      qs <- vapply(probs, function(pr) {
        at <- which.max(cw >= pr * total)
        if (cw[at] < pr * total) NA_real_ else yd[at]
      }, 0.0)
      pickands_from_quantiles(qs[1], qs[2], qs[3])
    }, 0.0)
    c(score = mean((local - marginal)^2), undefined = mean(!is.finite(local)))
  }, numeric(2))
  list(scores = out["score", ], undefined = out["undefined", ])
}

## --------------------------------------------------------------- (2) qSIS ---
qa_sis_scores_multi <- function(x, y, taus, df = 3) {
  n <- nrow(x); p <- ncol(x)
  qm <- vapply(taus, function(tt) as.numeric(quantile(y, tt, type = 7)), 0.0)
  out <- matrix(NA_real_, p, length(taus))
  for (j in seq_len(p)) {
    B <- cbind(1, bs(x[, j], df = df))   # built once, reused for every tau
    for (t in seq_along(taus)) {
      fit <- tryCatch(rq.fit.br(B, y, tau = taus[t]), error = function(e) NULL)
      if (is.null(fit)) next
      out[j, t] <- mean((B %*% fit$coefficients - qm[t])^2)
    }
  }
  out
}

## ------------------------------------------------------------------ check ---
set.seed(20260814)
z <- ar1_gaussian(n, p, 0.25)
u <- rank_uniform_matrix(z)
y <- exp(rnorm(n) * 2) + 1
k <- floor(0.072 * n)

ref_yu  <- yu_score_matrix(u, y, k, h = 1)
new_yu  <- yu_score_matrix_fast(u, y, k, h = 1)
cat("YU   identical scores:", identical(ref_yu$scores, new_yu$scores),
    " max abs diff:", max(abs(ref_yu$scores - new_yu$scores), na.rm = TRUE),
    " same ranking:", identical(order(-ref_yu$scores, seq_len(p)),
                                order(-new_yu$scores, seq_len(p))), "\n")

ref_q <- vapply(TAUS, function(tt) qa_sis_scores(u, y, tau = tt), numeric(p))
new_q <- qa_sis_scores_multi(u, y, TAUS)
cat("qSIS identical scores:", identical(ref_q, new_q),
    " max abs diff:", max(abs(ref_q - new_q), na.rm = TRUE), "\n")

## ------------------------------------------------------------------ timing --
tt <- function(e) unname(system.time(e)[3])
t_yu_ref <- tt(yu_score_matrix(u, y, k, h = 1))
t_yu_new <- tt(yu_score_matrix_fast(u, y, k, h = 1))
t_q_ref  <- tt(for (x in TAUS) qa_sis_scores(u, y, tau = x))
t_q_new  <- tt(qa_sis_scores_multi(u, y, TAUS))
cat(sprintf("n=%d p=%d\n  YU   %6.2f s -> %6.2f s  (%.2fx)\n  qSIS %6.2f s -> %6.2f s  (%.2fx)\n  total des deux %6.2f s -> %6.2f s (%.2fx)\n",
            n, p, t_yu_ref, t_yu_new, t_yu_ref / t_yu_new,
            t_q_ref, t_q_new, t_q_ref / t_q_new,
            t_yu_ref + t_q_ref, t_yu_new + t_q_new,
            (t_yu_ref + t_q_ref) / (t_yu_new + t_q_new)))
