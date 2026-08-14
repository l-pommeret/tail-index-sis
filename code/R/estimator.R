# Primary local-fraction Hill estimator following Gardes--Podgorny Definition 5.

pit_rank <- function(x) rank(x, ties.method = "average") / (length(x) + 1)

local_hill_r <- function(u, y, eval_u, h, alpha) {
  stopifnot(length(u) == length(y), h > 0, h < 1, alpha > 0, alpha < 1,
            all(is.finite(u)), all(is.finite(y)), all(y > 0))
  ans <- numeric(length(eval_u))
  local_n <- integer(length(eval_u))
  local_k <- integer(length(eval_u))
  underpopulated <- logical(length(eval_u))
  for (g in seq_along(eval_u)) {
    yy <- y[abs(u - eval_u[g]) <= h]
    m <- length(yy)
    k <- as.integer(floor(alpha * m))
    local_n[g] <- m
    local_k[g] <- k
    if (alpha * m <= 1 || k < 1L || k >= m) {
      ans[g] <- NA_real_
      underpopulated[g] <- TRUE
    } else {
      z <- sort(log(yy), decreasing = TRUE)
      ans[g] <- mean(z[seq_len(k)] - z[k + 1L])
    }
  }
  list(estimate = ans, local_n = local_n, local_k = local_k,
       underpopulated = underpopulated)
}

score_coordinate_r <- function(x, y, h, alpha, epsilon = 0.05) {
  u <- pit_rank(x)
  keep <- u >= epsilon & u <= 1 - epsilon
  fit <- local_hill_r(u, y, u[keep], h, alpha)
  valid <- is.finite(fit$estimate)
  c(score = if (any(valid)) mean(fit$estimate[valid]) else Inf,
    n_eval = sum(valid),
    under_rate = mean(fit$underpopulated), mean_local_n = mean(fit$local_n),
    mean_local_k = mean(fit$local_k))
}

score_matrix_r <- function(x, y, h, alpha, epsilon = 0.05) {
  stopifnot(is.matrix(x), nrow(x) == length(y))
  out <- vapply(seq_len(ncol(x)), function(j)
    score_coordinate_r(x[, j], y, h, alpha, epsilon), numeric(5))
  t(out)
}
