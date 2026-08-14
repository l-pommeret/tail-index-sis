ar1_gaussian <- function(n, p, rho, block_size = 256L) {
  stopifnot(n >= 2, p >= 1, abs(rho) < 1, block_size >= 1)
  z <- matrix(NA_real_, n, p)
  z[, 1L] <- rnorm(n)
  if (p > 1L) {
    innovation_sd <- sqrt(1 - rho^2)
    for (start in seq.int(2L, p, by = block_size)) {
      end <- min(p, start + block_size - 1L)
      for (j in start:end)
        z[, j] <- rho * z[, j - 1L] + innovation_sd * rnorm(n)
    }
  }
  z
}

rank_uniform_matrix <- function(z) {
  n <- nrow(z)
  apply(z, 2L, rank, ties.method = "average") / (n + 1)
}

population_uniform_matrix <- function(z) {
  stopifnot(is.matrix(z), all(is.finite(z)))
  pnorm(z)
}

gamma_model <- function(u, model = "A", signal = 1) {
  stopifnot(is.matrix(u), ncol(u) >= 4)
  if (model %in% c("A", "B", "F"))
    return(0.5 * exp(-signal * rowSums(u[, 1:4, drop = FALSE])))
  if (model == "C")
    return(0.18 + 0.32 * exp(-8 * rowSums((u[, 1:4, drop = FALSE] -
                                           matrix(c(.2, .4, .6, .8), nrow(u), 4,
                                                  byrow = TRUE))^2)))
  if (model == "D") {
    # Strictly decreasing in every active coordinate, with genuine interactions.
    return(0.55 * exp(-0.35 * (u[,1] + u[,2] + u[,3] + u[,4]) -
                        0.20 * u[,1] * u[,2] - 0.15 * u[,3] * u[,4]))
  }
  if (model == "E") return(0.35 - 0.20 * (u[,1] - u[,2])^2)
  if (model == "N") return(rep(0.3, nrow(u)))
  stop("unknown model")
}

wang_tsai_response <- function(gamma, m = 0.25, inactive_u = NULL,
                               model = "A") {
  stopifnot(all(is.finite(gamma)), all(gamma > 0), m >= 0)
  v <- runif(length(gamma))
  log_y <- gamma * log((m + 1 - m * v) / v)
  # Model F perturbs only the second-order/body component, not the tail index.
  if (model == "F") {
    stopifnot(!is.null(inactive_u))
    # Positive multiplicative perturbation on log(Y), asymptotic to an
    # inactive-covariate-dependent scale factor on Y. It preserves gamma.
    log_y <- log_y * (1 + 0.20 * (inactive_u - 0.5) / (1 + log_y))
  }
  y <- exp(log_y)
  if (any(!is.finite(y)) || any(y <= 1)) stop("Wang--Tsai invariant Y > 1 failed")
  y
}

simulate_dataset <- function(n, p, rho, model, signal = 1, seed,
                             block_size = 256L) {
  set.seed(seed)
  z <- ar1_gaussian(n, p, rho, block_size)
  # The data-generating truth is a function of the population PIT.  The
  # estimator deliberately does not receive `u`: it computes empirical ranks
  # from `x`, as required by the rank-based procedure.
  u <- population_uniform_matrix(z)
  gamma <- gamma_model(u, model, signal)
  inactive_u <- if (p >= 5) u[,5] else u[,p]
  y <- wang_tsai_response(gamma, inactive_u = inactive_u, model = model)
  list(x = z, u = u, y = y, gamma = gamma)
}
