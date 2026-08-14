## Quantile-adaptive model-free screening of He, Wang and Hong (2013, AoS).
## For each coordinate j, fit the tau-th conditional quantile of Y on a cubic
## B-spline basis of X_j and rank coordinates by the marginal utility
##   w_j = n^{-1} sum_i { fhat_j(X_ij) - Fy^{-1}(tau) }^2,
## descending (large utility ranks first).  Implementation follows their
## Section 2.2 with df = 3 spline basis, as in their simulations.

suppressPackageStartupMessages({
  library(quantreg)
  library(splines)
})

qa_sis_scores <- function(x, y, tau = 0.95, df = 3) {
  n <- nrow(x)
  q_marginal <- as.numeric(quantile(y, tau, type = 7))
  vapply(seq_len(ncol(x)), function(j) {
    B <- bs(x[, j], df = df)
    fit <- tryCatch(rq.fit.br(cbind(1, B), y, tau = tau),
                    error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    fhat <- cbind(1, B) %*% fit$coefficients
    mean((fhat - q_marginal)^2)
  }, numeric(1))
}
