## Simulation models, Draft 5.
##
## Family A -- tail-index structure.  A1, A2, A3 share the Gaussian AR(1)
## covariate design, the slowly varying factor ell_1 of Gardes-Podgorny, and
## the tail-index active set A_gamma = {1,2,3,4}; they differ only in the
## tail-index function:
##   A1  gamma = 0.50 exp(-s4)                             (reference)
##   A2  gamma = 0.50 exp(-0.5 s4)                         (weak signal)
##   A3  gamma = 0.50 exp(-0.80 s4 - 0.20 u1 u2 - 0.15 u3 u4)  (interactions)
## with s4 = u1 + u2 + u3 + u4, u = Phi(z).
##
## Family B -- finite-quantile scale contamination.  B1 has the same tail
## index as A1.  A latent standard normal factor F, independent of the active
## coordinates, enters the response through
##   log ell_2 = -1/(2s) + 0.5 (Phi(F) - 0.5),
## so the persistent scale factor lies in [e^-0.25, e^0.25].  F is not
## observed; q_scale = 20 proxies occupy coordinates A_SCALE = {5,...,24}:
##   Z_j = 0.7 F + sqrt(1 - 0.49) eps_j,   j in A_SCALE,
## each marginally N(0,1), pairwise intra-block correlation 0.49.  The B1
## covariate vector is the base AR(1) vector with coordinates 5..24 replaced
## by the proxy block; coordinates 1..4 and 25..p keep the base AR(1) law,
## and the AR(1) chain beyond 24 continues from the BASE coordinate 24, not
## from the overwritten proxy (full and streaming generators agree).

source("code/R/generate.R")   # ar1_gaussian, population_uniform_matrix

B1_LAMBDA     <- 0.7
B1_SCALE_SIZE <- 20L
A_GAMMA <- 1:4
A_SCALE <- 5:(4L + B1_SCALE_SIZE)

gamma_model <- function(u, model) {
  s4 <- rowSums(u[, 1:4, drop = FALSE])
  switch(model,
    A1 = 0.50 * exp(-s4),
    A2 = 0.50 * exp(-0.5 * s4),
    A3 = 0.50 * exp(-0.80 * s4 - 0.20 * u[, 1] * u[, 2] -
                    0.15 * u[, 3] * u[, 4]),
    B1 = 0.50 * exp(-s4),
    stop("unknown model"))
}

sv_l1 <- function(u) 2.5 * rowSums(u[, 5:8, drop = FALSE])

## B1 proxy block: overwrite coordinates A_SCALE with latent-factor proxies.
b1_scale_block <- function(z, f, lambda = B1_LAMBDA) {
  n <- nrow(z)
  for (j in A_SCALE)
    z[, j] <- lambda * f + sqrt(1 - lambda^2) * rnorm(n)
  z
}

response5 <- function(gamma, u, model, V, f = NULL) {
  stopifnot(all(gamma > 0))
  ell <- if (model == "B1") {
    stopifnot(!is.null(f))
    exp(-V / 2) * exp(0.5 * (pnorm(f) - 0.5))
  } else {
    1 / (1 + exp(sv_l1(u) - pmin(1 / V, 700)))
  }
  y <- V^(-gamma) * ell
  stopifnot(all(is.finite(y)), all(y > 0))
  y
}

simulate_dataset5 <- function(n, p, rho, model, seed, block_size = 256L) {
  stopifnot(p >= max(A_SCALE))
  set.seed(seed)
  z <- ar1_gaussian(n, p, rho, block_size)
  f <- NULL
  if (model == "B1") {
    f <- rnorm(n)
    z <- b1_scale_block(z, f)
  }
  u <- population_uniform_matrix(z)
  gamma <- gamma_model(u, model)
  y <- response5(gamma, u, model, V = runif(n), f = f)
  list(z = z, u = u, y = y, gamma = gamma, f = f)
}

## Streaming version for tuning sweeps.  The initial block covers all
## model-relevant coordinates (max(A_SCALE)).  For B1 the base AR(1) state at
## coordinate max(A_SCALE) is preserved BEFORE the proxy block overwrites
## coordinates 5..24, so that coordinates beyond the block continue the base
## AR(1) chain exactly as in simulate_dataset5.
simulate_scores_streaming5 <- function(n, p, rho, model, seed, settings,
                                       epsilon = .05, block_size = 256L) {
  ## settings: data.frame with columns h, alpha; one score vector per row.
  n0 <- max(A_SCALE)
  stopifnot(p >= n0)
  set.seed(seed)
  innovation_sd <- sqrt(1 - rho^2)
  S <- matrix(NA_real_, p, nrow(settings))
  initial <- matrix(NA_real_, n, n0)
  initial[, 1L] <- rnorm(n)
  for (j in 2:n0)
    initial[, j] <- rho * initial[, j - 1L] + innovation_sd * rnorm(n)
  previous_ar <- initial[, n0]          # base AR(1) state at coordinate n0
  f <- NULL
  if (model == "B1") {
    f <- rnorm(n)
    initial <- b1_scale_block(initial, f)
  }
  u_initial <- population_uniform_matrix(initial)
  gamma <- gamma_model(u_initial, model)
  y <- response5(gamma, u_initial, model, V = runif(n), f = f)
  score_block <- function(z, indices) {
    for (k in seq_len(nrow(settings))) {
      dg <- vapply(seq_along(indices), function(q)
        score_coordinate_cpp(z[, q], y, settings$h[k], settings$alpha[k],
                             epsilon), numeric(5))
      S[indices, k] <<- dg["score", ]
    }
  }
  score_block(initial, seq_len(n0))
  previous <- previous_ar               # continue from the base chain
  rm(initial, u_initial)
  if (p > n0) {
    for (start in seq.int(n0 + 1L, p, by = block_size)) {
      end <- min(p, start + block_size - 1L)
      z <- matrix(NA_real_, n, end - start + 1L)
      for (q in seq_len(ncol(z))) {
        z[, q] <- rho * previous + innovation_sd * rnorm(n)
        previous <- z[, q]
      }
      score_block(z, start:end)
    }
  }
  S
}
