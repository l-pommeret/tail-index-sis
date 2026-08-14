## Suite de modèles figée (Draft 4).
##
## Écrite à côté de generate3.R et non à sa place : generate3.R reste la
## référence Draft 3, pour que reproduce_draft3.sh continue de reproduire les
## résultats publiés et que results/draft3/SHA256SUMS.txt reste vérifiable.
##
## Différences avec Draft 3, décidées d'après les mesures consignées dans
## results/decisions_modeles.md :
##
##   M1  inchangé.   gamma = 0.50 exp(-s4)
##   M2  facteur d'échelle en AMAS CORRÉLÉ (D4). gamma inchangé ; A_scale passe
##       de {5..8} à {5..24}, les 20 coordonnées étant des proxys d'un facteur
##       latent unique F de charge lambda = 0.7, et log ell = kappa_c (Phi(F)-1/2).
##   M3  inchangé.   gamma = 0.50 exp(-0.5 s4)
##   M4  taux de décroissance porté de 0.35 à 0.80 (D1) ; interactions,
##       constante de tête et facteur à variation lente inchangés.
##
## Pourquoi l'amas pour M2.  Le facteur publié log ell = kappa sum_{A_scale}
## (u_j - 1/2) lie l'attractivité d'un leurre (proportionnelle à kappa) au dégât
## infligé au Hill local (proportionnel à kappa^2 s) : le nombre de covariables
## d'échelle ne peut donc pas servir de cadran.  Avec un facteur latent unique,
## Var(log ell) = kappa_c^2/12 quels que soient s et lambda, donc la difficulté
## est figée pendant que chaque proxy garde une corrélation marginale lambda
## avec la nuisance.
##
## Pourquoi kappa_c = 0.4772.  sd(log ell) = kappa_c/sqrt(12) EXACTEMENT, sans
## dépendance en rho, contrairement à la version publiée dont la somme est
## gonflée par la corrélation AR(1) du bloc.  La valeur retenue est calée sur la
## sd mesurée du M2 publié au point de référence de la Draft 3 (rho = 0.25),
## soit 0.1377, et gardée FIXE ensuite : la magnitude de la nuisance ne varie
## donc plus avec rho, ce qui supprime un facteur de confusion que la version
## publiée introduisait.

source("code/R/generate.R")   # ar1_gaussian, population_uniform_matrix

## ------------------------------------------------------------------ M2 ------
M2_LAMBDA  <- as.numeric(Sys.getenv("M2_LAMBDA", "0.7"))
M2_S       <- as.integer(Sys.getenv("M2_S", "20"))
M2_KAPPA_C <- as.numeric(Sys.getenv("M2_KAPPA_C", "0.4772"))
A_GAMMA <- 1:4                       # actives en indice de queue
A_SCALE <- 5:(4L + M2_S)             # actives en échelle (M2 uniquement)

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

## Facteur à variation lente de Gardes-Podgorny, pour M1, M3, M4 : il tend vers
## 1 quand V -> 0, donc son effet s'évanouit dans la queue.
sv_l1_4 <- function(u) 2.5 * rowSums(u[, 5:8, drop = FALSE])

## Tirage du bloc d'échelle de M2 : z_j = lambda F + sqrt(1-lambda^2) eps_j.
## Les colonnes sont écrasées après le tirage AR(1), donc la corrélation du bloc
## avec ses deux colonnes voisines n'est par construction pas rho ; tout ce qui
## est hors du bloc est intact.
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
## initial doit couvrir toutes les coordonnées qui entrent dans le modèle, donc
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
