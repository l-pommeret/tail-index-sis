## Suite de modèles de simulation (Draft 4).
##
## Écrite à côté de generate3.R et non à sa place : generate3.R reste la
## référence Draft 3, pour que reproduce_draft3.sh continue de reproduire les
## résultats publiés et que results/draft3/SHA256SUMS.txt reste vérifiable.
##
## Ensemble actif en indice de queue : A_gamma = {1,2,3,4} dans les quatre
## modèles.  Ensemble actif en échelle : A_scale = {5,...,24}, pour M2
## uniquement ; ces coordonnées déplacent les quantiles conditionnels finis sans
## toucher l'indice de queue.
##
##   M1  gamma = 0.50 exp(-s4)
##   M2  gamma = 0.50 exp(-s4), avec facteur d'échelle
##   M3  gamma = 0.50 exp(-0.5 s4)
##   M4  gamma = 0.55 exp(-0.80 s4 - 0.20 u1 u2 - 0.15 u3 u4)
##
## où s4 = u1 + u2 + u3 + u4 et u = Phi(z) est la transformation intégrale de
## probabilité de population.
##
## Facteurs à variation lente.  M1, M3 et M4 portent le facteur de
## Gardes-Podgorny ell = 1/(1 + exp(l1(u) - 1/V)), qui tend vers 1 quand V -> 0
## et dont l'effet s'évanouit donc dans la queue.  M2 porte un facteur d'échelle
## qui persiste : un facteur latent unique F, observé à travers les s = 20
## coordonnées de A_scale, chacune de charge lambda,
##
##   z_j    = lambda F + sqrt(1 - lambda^2) eps_j,     j dans A_scale
##   log ell = kappa_c (Phi(F) - 1/2)      (plus le terme -V/2)
##
## Var(log ell) = kappa_c^2 / 12 exactement, sans dépendance en s, en lambda ni
## en rho : la magnitude de la nuisance est donc constante sur tout l'espace de
## design, et le nombre de coordonnées d'échelle ainsi que leur qualité de proxy
## se règlent indépendamment d'elle.

source("code/R/generate.R")   # ar1_gaussian, population_uniform_matrix

M2_LAMBDA  <- as.numeric(Sys.getenv("M2_LAMBDA", "0.7"))
M2_S       <- as.integer(Sys.getenv("M2_S", "20"))
M2_KAPPA_C <- as.numeric(Sys.getenv("M2_KAPPA_C", "0.4772"))
A_GAMMA <- 1:4
A_SCALE <- 5:(4L + M2_S)

gamma_model4 <- function(u, model) {
  s4 <- rowSums(u[, 1:4, drop = FALSE])
  switch(model,
    M1 = 0.50 * exp(-s4),
    M2 = 0.50 * exp(-s4),
    M3 = 0.50 * exp(-0.5 * s4),
    M4 = 0.55 * exp(-0.80 * s4 - 0.20 * u[, 1] * u[, 2] -
                    0.15 * u[, 3] * u[, 4]),
    stop("unknown model"))
}

sv_l1_4 <- function(u) 2.5 * rowSums(u[, 5:8, drop = FALSE])

## Bloc d'échelle de M2.  Les colonnes sont écrasées après le tirage AR(1), donc
## la corrélation du bloc avec ses deux colonnes voisines n'est par construction
## pas rho ; tout ce qui est hors du bloc est intact.
m2_scale_block <- function(z, f, lambda = M2_LAMBDA) {
  n <- nrow(z)
  for (j in A_SCALE)
    z[, j] <- lambda * f + sqrt(1 - lambda^2) * rnorm(n)
  z
}

response4 <- function(gamma, u, model, V, f = NULL) {
  stopifnot(all(gamma > 0))
  ell <- if (model == "M2") {
    stopifnot(!is.null(f))
    exp(-V / 2) * exp(M2_KAPPA_C * (pnorm(f) - 0.5))
  } else {
    1 / (1 + exp(sv_l1_4(u) - pmin(1 / V, 700)))
  }
  y <- V^(-gamma) * ell
  stopifnot(all(is.finite(y)), all(y > 0))
  y
}

simulate_dataset4 <- function(n, p, rho, model, seed, block_size = 256L) {
  stopifnot(p >= max(A_SCALE))
  set.seed(seed)
  z <- ar1_gaussian(n, p, rho, block_size)
  f <- NULL
  if (model == "M2") {
    f <- rnorm(n)
    z <- m2_scale_block(z, f)
  }
  u <- population_uniform_matrix(z)
  gamma <- gamma_model4(u, model)
  y <- response4(gamma, u, model, V = runif(n), f = f)
  list(z = z, u = u, y = y, gamma = gamma, f = f)
}

## Version en flux, pour les balayages de réglage en grande dimension.  Le bloc
## initial couvre toutes les coordonnées qui entrent dans le modèle, donc
## max(A_SCALE) et non 8 comme dans generate3.R.
simulate_score_streaming4 <- function(n, p, rho, model, seed, h, alpha,
                                      epsilon = .05, block_size = 256L) {
  n0 <- max(A_SCALE)
  stopifnot(p >= n0)
  set.seed(seed)
  innovation_sd <- sqrt(1 - rho^2)
  scores <- numeric(p)
  initial <- matrix(NA_real_, n, n0)
  initial[, 1L] <- rnorm(n)
  for (j in 2:n0)
    initial[, j] <- rho * initial[, j - 1L] + innovation_sd * rnorm(n)
  f <- NULL
  if (model == "M2") {
    f <- rnorm(n)
    initial <- m2_scale_block(initial, f)
  }
  u_initial <- population_uniform_matrix(initial)
  gamma <- gamma_model4(u_initial, model)
  y <- response4(gamma, u_initial, model, V = runif(n), f = f)
  score_block <- function(z, indices) {
    dg <- vapply(seq_along(indices), function(q)
      score_coordinate_cpp(z[, q], y, h, alpha, epsilon), numeric(5))
    scores[indices] <<- dg["score", ]
  }
  score_block(initial, seq_len(n0))
  previous <- initial[, n0]
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
  scores
}
