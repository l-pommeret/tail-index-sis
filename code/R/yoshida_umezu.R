# Yoshida--Umezu conditional-Pickands screening statistic (their eqs. 2--4).
# The implementation sorts Y once and evaluates weighted empirical survival
# quantiles. Larger utilities rank first, unlike the proposed local-Hill score.

weighted_upper_quantile <- function(y_desc, w_desc, tail_probability) {
  total <- sum(w_desc)
  if (!is.finite(total) || total <= 0) return(NA_real_)
  at <- which(cumsum(w_desc) >= tail_probability * total)[1L]
  if (is.na(at)) return(NA_real_)
  y_desc[at]
}

pickands_from_quantiles <- function(q1, q2, q4) {
  ratio <- (q1 - q2) / (q2 - q4)
  if (!is.finite(ratio) || ratio <= 0) return(NA_real_)
  log(ratio) / log(2)
}

yu_score_coordinate <- function(u, y, k, h = 1, eval_z = seq(.02, .98, length.out=25)) {
  n <- length(y)
  stopifnot(length(u) == n, k >= 1L, 4L*k < n, h > 0)
  oy <- order(y, decreasing=TRUE)
  yd <- y[oy]; ud <- u[oy]
  probs <- c(k, 2*k, 4*k) / n
  marginal <- pickands_from_quantiles(yd[k], yd[2*k], yd[4*k])
  local <- vapply(eval_z, function(z) {
    zz <- (ud-z)/h
    w <- (1-zz^2) * (abs(zz) <= 1) # Epanechnikov, normalization cancels
    qs <- vapply(probs, function(pr) weighted_upper_quantile(yd,w,pr), 0.0)
    pickands_from_quantiles(qs[1],qs[2],qs[3])
  }, 0.0)
  c(score=mean((local-marginal)^2), undefined=mean(!is.finite(local)))
}

yu_score_matrix <- function(u, y, k, h=1, eval_z=seq(.02,.98,length.out=25)) {
  out <- vapply(seq_len(ncol(u)), function(j)
    yu_score_coordinate(u[,j],y,k,h,eval_z), numeric(2))
  list(scores=out["score",], undefined=out["undefined",])
}
