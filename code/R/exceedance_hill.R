## Exceedance-regression screen: standardise once, then average.
##
## For a conditional Pareto-type response, log(Y / Q(tau0|u)) given an
## exceedance has conditional mean gamma(u).  So:
##   1. fit ONE intermediate conditional quantile Q(tau0 | X_j) by quantile
##      regression on a cubic B-spline basis -- dividing by it removes any
##      covariate effect that acts on the scale rather than on the tail;
##   2. among the exceedances only, regress Z = log(Y/Qhat) on the same basis by
##      least squares; the fitted curve estimates gamma_j(.), and every fitted
##      value averages of order (1-tau0)n log-spacings.
##
## The contrast with code/R/qr_slope.R is the point: that screen cancelled the
## scale by DIFFERENCING two quantile fits, which amplifies variance and made it
## the worst screen tested; this one cancels the scale by DIVISION and then
## AVERAGES, which is where the Hill estimator's efficiency comes from.
##
## Cost is one rq fit per coordinate instead of the four quantile SIS uses, plus
## a least-squares fit on the exceedances, so it is cheaper than quantile SIS.
##
## Three readings are returned, since level-versus-variation is not settled:
##   g_mean  average of ghat along the fibre, ranked ASCENDING
##   g_sd    dispersion of ghat along the fibre, ranked DESCENDING
##   g_F     F statistic for flatness of the spline terms, ranked DESCENDING;
##           under an inactive coordinate the envelope is flat, so this one has
##           a known null distribution and yields scores comparable across
##           coordinates.

suppressPackageStartupMessages({
  library(quantreg)
  library(splines)
})

exceedance_hill_scores <- function(x, y, tau0 = 0.90, df = 3) {
  n <- nrow(x); p <- ncol(x)
  ly <- log(y)
  out <- matrix(NA_real_, p, 3L,
                dimnames = list(NULL, c("g_mean", "g_sd", "g_F")))
  for (j in seq_len(p)) {
    B <- cbind(1, bs(x[, j], df = df))
    fit <- tryCatch(rq.fit.br(B, ly, tau = tau0), error = function(e) NULL)
    if (is.null(fit)) next
    lq <- as.numeric(B %*% fit$coefficients)     # fitted log conditional quantile
    z <- ly - lq
    ex <- which(z > 0)                            # exceedances
    if (length(ex) < 5 * (df + 1)) next
    Be <- B[ex, , drop = FALSE]
    ze <- z[ex]
    ls <- tryCatch(lm.fit(Be, ze), error = function(e) NULL)
    if (is.null(ls) || anyNA(ls$coefficients)) next
    ghat <- as.numeric(B %*% ls$coefficients)     # gamma_j(.) over the fibre
    ## F test of flatness: spline terms jointly zero
    rss1 <- sum(ls$residuals^2)
    rss0 <- sum((ze - mean(ze))^2)
    dfr <- length(ex) - (df + 1L)
    out[j, "g_mean"] <- mean(ghat)
    out[j, "g_sd"] <- sd(ghat)
    out[j, "g_F"] <- if (rss1 > 0 && dfr > 0)
      ((rss0 - rss1) / df) / (rss1 / dfr) else NA_real_
  }
  out
}
