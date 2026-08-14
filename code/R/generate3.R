## Draft-3 simulation designs.
## Identical to Draft 2 (generate2.R) except for the slowly varying factor of
## Model M2, recalibrated as
##   l2(s, u) = exp{-1/(2s)} * v_kappa(u),
##   v_kappa(u) = exp{ kappa * sum_{j=5}^{8} (u_j - 1/2) },
## so that e^{-2 kappa} <= v_kappa <= e^{2 kappa}: strictly positive,
## uniformly bounded, symmetric in the inactive variables U5..U8, and
## controlled by the single strength parameter kappa.  The tail-index active
## set is A_gamma = {1,2,3,4}; A_scale = {5,6,7,8} moves finite conditional
## quantiles only.  kappa is fixed by the pre-specified pilot protocol
## (KAPPA env var; pilot seed stream disjoint from all final streams).

source("code/R/generate.R")  # ar1_gaussian, population_uniform_matrix

KAPPA <- as.numeric(Sys.getenv("KAPPA", "0.35"))

gamma_model3 <- function(u, model) {
  s4 <- rowSums(u[, 1:4, drop = FALSE])
  switch(model,
    M1 = 0.5 * exp(-s4),
    M2 = 0.5 * exp(-s4),
    M3 = 0.5 * exp(-0.5 * s4),
    M4 = 0.55 * exp(-0.35 * s4 - 0.20 * u[, 1] * u[, 2] -
                    0.15 * u[, 3] * u[, 4]),
    stop("unknown model"))
}

## l1 keeps the Gardes-Podgorny form with v(u) = 2.5 * sum(u5..u8) in (0,10).
sv_l1 <- function(u) 2.5 * rowSums(u[, 5:8, drop = FALSE])
## l2 scale factor, Draft-3 recalibration.
sv_l2 <- function(u, kappa = KAPPA)
  exp(kappa * rowSums(u[, 5:8, drop = FALSE] - 0.5))

response3 <- function(gamma, u, model, V = runif(length(gamma))) {
  stopifnot(all(gamma > 0))
  s <- 1 / V
  ell <- if (model == "M2") exp(-V / 2) * sv_l2(u)
         else 1 / (1 + exp(sv_l1(u) - pmin(s, 700)))
  y <- V^(-gamma) * ell
  stopifnot(all(is.finite(y)), all(y > 0))
  y
}

simulate_dataset3 <- function(n, p, rho, model, seed, block_size = 256L) {
  set.seed(seed)
  z <- ar1_gaussian(n, p, rho, block_size)
  u <- population_uniform_matrix(z)
  gamma <- gamma_model3(u, model)
  y <- response3(gamma, u, model)
  list(z = z, u = u, y = y, gamma = gamma)
}

simulate_score_streaming3 <- function(n, p, rho, model, seed, h, alpha,
                                      epsilon = .05, block_size = 256L) {
  stopifnot(p >= 8L)
  set.seed(seed)
  innovation_sd <- sqrt(1 - rho^2)
  scores <- numeric(p)
  initial <- matrix(NA_real_, n, 8L)
  initial[, 1L] <- rnorm(n)
  for (j in 2:8)
    initial[, j] <- rho * initial[, j - 1L] + innovation_sd * rnorm(n)
  u_initial <- population_uniform_matrix(initial)
  gamma <- gamma_model3(u_initial, model)
  y <- response3(gamma, u_initial, model)
  score_block <- function(z, indices) {
    dg <- vapply(seq_along(indices), function(q)
      score_coordinate_cpp(z[, q], y, h, alpha, epsilon), numeric(5))
    scores[indices] <<- dg["score", ]
  }
  score_block(initial, 1:8)
  previous <- initial[, 8L]
  rm(initial, u_initial)
  if (p > 8L) {
    for (start in seq.int(9L, p, by = block_size)) {
      end <- min(p, start + block_size - 1L)
      z <- matrix(NA_real_, n, end - start + 1L)
      for (q in seq_len(ncol(z))) {
        z[, q] <- rho * previous + innovation_sd * rnorm(n)
        previous <- z[, q]
      }
      score_block(z, start:end)
    }
  }
  scores
}
