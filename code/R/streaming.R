source("code/R/generate.R")

simulate_score_streaming <- function(n, p, rho, model, signal, seed, h, alpha,
                                     epsilon = .05, block_size = 256L) {
  stopifnot(p >= 5L, block_size >= 1L)
  set.seed(seed)
  innovation_sd <- sqrt(1 - rho^2)
  scores <- numeric(p)
  under_rate <- numeric(p)
  mean_local_n <- numeric(p)
  mean_local_k <- numeric(p)

  # The response requires only active columns 1:4 and column 5 for model F.
  # Always generate these five first, independently of the processing block size.
  initial <- matrix(NA_real_, n, 5L)
  initial[, 1L] <- rnorm(n)
  for (j in 2:5)
    initial[, j] <- rho * initial[, j - 1L] + innovation_sd * rnorm(n)
  # Construct the response from population PITs.  score_coordinate_cpp below
  # still forms empirical ranks separately for every observed covariate.
  u_initial <- population_uniform_matrix(initial)
  gamma <- gamma_model(u_initial, model, signal)
  y <- wang_tsai_response(gamma, inactive_u = u_initial[,5L], model = model)

  score_block <- function(z, indices) {
    dg <- vapply(seq_along(indices), function(q)
      score_coordinate_cpp(z[, q], y, h, alpha, epsilon), numeric(5))
    scores[indices] <<- dg["score",]
    under_rate[indices] <<- dg["under_rate",]
    mean_local_n[indices] <<- dg["mean_local_n",]
    mean_local_k[indices] <<- dg["mean_local_k",]
  }
  score_block(initial, 1:5)
  previous <- initial[,5L]
  rm(initial, u_initial)

  if (p > 5L) {
    for (start in seq.int(6L, p, by = block_size)) {
      end <- min(p, start + block_size - 1L)
      z <- matrix(NA_real_, n, end - start + 1L)
      for (q in seq_len(ncol(z))) {
        z[,q] <- rho * previous + innovation_sd * rnorm(n)
        previous <- z[,q]
      }
      score_block(z, start:end)
    }
  }
  list(scores=scores, under_rate=under_rate, mean_local_n=mean_local_n,
       mean_local_k=mean_local_k, y_min=min(y), y_finite=all(is.finite(y)),
       gamma_range=range(gamma))
}
