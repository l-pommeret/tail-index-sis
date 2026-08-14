## Extremal quantile-regression slope screen.
##
## For a heavy-tailed response, the conditional quantile satisfies
##   log Q(tau | X_j = u) = gamma_j(u) * log{1/(1-tau)} + log C_j(u) + o(1),
## a straight line in s = log{1/(1-tau)} whose SLOPE is the conditional tail
## index and whose INTERCEPT is the scale.  Quantile SIS reads the line at a
## single tau, so it cannot tell the two apart and a scale-only covariate fools
## it; the local-Hill screen targets the slope but estimates it from about two
## independent rank windows.  This screen fits the same quantile-regression
## machinery quantile SIS uses, at L extreme levels, and keeps only the slope:
## a covariate that multiplies Y by a factor free of tau shifts the intercept
## and leaves the slope exactly unchanged, so it is cancelled by construction.
##
## Conditioning on X_j = u alone leaves a mixture over the other coordinates,
## whose tail index is the maximum along the fibre, so gamma_j(.) estimates the
## same upper envelope Psi_j(.) the local-Hill score targets: same estimand,
## different estimator.
##
## Two rankings are returned, as the level/variation contrast is not settled:
##   slope_mean  average slope along the fibre, ranked ASCENDING (flat at the
##               global maximum for an inactive coordinate, depressed for an
##               active one) -- the analogue of the current score;
##   slope_sd    dispersion of the slope along the fibre, ranked DESCENDING
##               (flat for an inactive coordinate) -- cancels any bias common
##               to all coordinates.

suppressPackageStartupMessages({
  library(quantreg)
  library(splines)
})

qr_slope_scores <- function(x, y, taus = c(0.90, 0.95, 0.975, 0.99), df = 3) {
  n <- nrow(x); p <- ncol(x)
  ly <- log(y)
  s <- log(1 / (1 - taus))
  sc <- s - mean(s)
  denom <- sum(sc^2)
  out <- matrix(NA_real_, p, 2L,
                dimnames = list(NULL, c("slope_mean", "slope_sd")))
  for (j in seq_len(p)) {
    B <- cbind(1, bs(x[, j], df = df))
    Q <- vapply(taus, function(tt) {
      fit <- tryCatch(rq.fit.br(B, ly, tau = tt), error = function(e) NULL)
      if (is.null(fit)) rep(NA_real_, n)
      else as.numeric(B %*% fit$coefficients)
    }, numeric(n))
    if (anyNA(Q)) next
    ## per-observation OLS slope of the L fitted log-quantiles on s
    slope <- as.numeric(Q %*% sc) / denom
    out[j, "slope_mean"] <- mean(slope)
    out[j, "slope_sd"] <- sd(slope)
  }
  out
}
