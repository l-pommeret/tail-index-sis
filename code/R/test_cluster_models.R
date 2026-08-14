## The M2 scale perturbation, in its correlated-cluster form, grafted onto the
## other models.
##
## M1, M3 and M4 carry the Gardes-Podgorny factor 1/(1+exp(l1(u) - 1/V)), which
## tends to 1 as V -> 0: the effect of coordinates 5..8 vanishes in the tail, so
## these models have no scale-active variables at all.  M2 alone carries
## exp(-V/2) * v_kappa(u), which persists.  Here the persistent factor is put on
## every model, in three forms:
##
##   publie   the model exactly as published
##   indep    s=4 scale covariates, log ell = kappa * sum(u_j - 1/2)   [= M2]
##   amas     s=20 proxies of one latent factor F, log ell = kappa_c (Phi(F)-1/2)
##            with kappa_c calibrated so Var(log ell) matches the indep s=4 arm
##
## gamma keeps each model's own form (code/R/generate3.R), so what varies across
## arms is only the nuisance, never the tail-index signal.
##
## usage: KAPPA=0.20 Rscript code/R/test_cluster_models.R OUT.rds MODELS [NREP] [CORES]
##   MODELS: comma separated, e.g. M3,M4

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R"); source("code/R/qa_sis.R")

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1L) args[1L] else "results/cluster_models/raw.rds"
MODELS <- if (length(args) >= 2L) strsplit(args[2L], ",")[[1]] else c("M3", "M4")
NREP  <- if (length(args) >= 3L) as.integer(args[3L]) else 40L
CORES <- if (length(args) >= 4L) as.integer(args[4L]) else 92L
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
KAP <- as.numeric(Sys.getenv("KAPPA", "0.20"))
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.05, 0.10, 0.15)
n <- 2000L; p <- 1000L; rho <- 0.25
DS <- c(4, 10, 20, 30, 50)
LAMBDAS <- c(0.7, 0.9)
S_CLUSTER <- 20L

draw <- function(seed, model, design, s, kappa, lambda = NA_real_) {
  set.seed(seed)
  z <- ar1_gaussian(n, p, rho)
  scale_idx <- 5:(4 + s)
  f <- NULL
  if (design == "amas") {
    f <- rnorm(n)
    for (j in scale_idx)
      z[, j] <- lambda * f + sqrt(1 - lambda^2) * rnorm(n)
  }
  u <- population_uniform_matrix(z)
  gamma <- gamma_model3(u, model)          # each model keeps its own gamma
  V <- runif(n)
  ## "publie" keeps the model's own slowly varying factor: l1 for M1/M3/M4,
  ## l2 for M2.  The other arms SUBSTITUTE l2 FOR l1, i.e. they use M2's
  ## branch of response3, ell = exp(-V/2) * v(u), with v built either from s
  ## independent scale coordinates or from the latent-factor cluster.
  log_v <- if (design == "amas") kappa * (pnorm(f) - 0.5)
           else kappa * rowSums(u[, scale_idx, drop = FALSE] - 0.5)
  ell <- if (design == "publie") {
    if (model == "M2") exp(-V / 2) * sv_l2(u)
    else 1 / (1 + exp(sv_l1(u) - pmin(1 / V, 700)))
  } else {
    exp(-V / 2) * exp(log_v)          # l2, substituted for l1
  }
  y <- V^(-gamma) * ell
  stopifnot(all(is.finite(y)), all(y > 0))
  list(z = z, y = y, scale_idx = scale_idx,
       sd_log_ell = if (design == "publie") NA_real_ else sd(log_v))
}

cal <- vapply(1:20, function(r)
  draw(61000019L + r * 401L, "M2", "indep", 4L, KAP)$sd_log_ell, numeric(1))
SD_TARGET <- mean(cal); KAPPA_C <- SD_TARGET * sqrt(12)
cat(sprintf("calibration : sd(log ell) reference = %.4f -> kappa_c = %.4f\n",
            SD_TARGET, KAPPA_C))

ARMS <- list()
for (m in MODELS) {
  ARMS[[length(ARMS) + 1L]] <- list(model = m, design = "publie", s = 4L,
    kappa = NA_real_, lambda = NA_real_, tag = sprintf("%s publie", m))
  ARMS[[length(ARMS) + 1L]] <- list(model = m, design = "indep", s = 4L,
    kappa = KAP, lambda = NA_real_,
    tag = sprintf("%s + echelle s=4 independante", m))
  for (l in LAMBDAS)
    ARMS[[length(ARMS) + 1L]] <- list(model = m, design = "amas",
      s = S_CLUSTER, kappa = KAPPA_C, lambda = l,
      tag = sprintf("%s + amas s=20, lambda=%.1f", m, l))
}

score_vec <- function(z, y, a, b) {
  h <- n^(-b) / 2; alpha <- n^(-a)
  vapply(seq_len(p), function(j)
    score_coordinate_cpp(z[, j], y, h, alpha, EPS), numeric(5))["score", ]
}

jobs <- do.call(rbind, lapply(seq_along(ARMS), function(i)
  data.frame(arm = i, r = seq_len(NREP))))

one <- function(ix) {
  A <- ARMS[[jobs$arm[ix]]]; r <- jobs$r[ix]
  seed <- 61000019L + jobs$arm[ix] * 100003L + r * 307L
  d <- draw(seed, A$model, A$design, A$s, A$kappa, A$lambda)
  base <- score_vec(d$z, d$y, ASTAR, BSTAR)
  R <- matrix(NA_real_, p, length(AGRID) * length(BGRID)); ii <- 0L
  for (a in AGRID) for (b in BGRID) {
    ii <- ii + 1L
    R[, ii] <- rank(score_vec(d$z, d$y, a, b), ties.method = "first",
                    na.last = TRUE)
  }
  agg <- apply(R, 1L, min)
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (n + 1))
  ords <- list(order(base, seq_len(p), na.last = TRUE),
               order(agg, seq_len(p)),
               order(-qa_sis_scores(uh, d$y, tau = 0.95), seq_len(p)))
  list(arm = jobs$arm[ix], replicate = r, seed = seed,
       sd_log_ell = d$sd_log_ell,
       rmax = vapply(ords, function(o) max(match(1:4, o)), numeric(1)),
       t4g = vapply(ords, function(o) sum(o[1:4] %in% 1:4), numeric(1)),
       t4s = vapply(ords, function(o) sum(o[1:4] %in% d$scale_idx), numeric(1)),
       t24s = vapply(ords, function(o) sum(o[1:24] %in% d$scale_idx), numeric(1)))
}

res <- parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = CORES,
                          mc.preschedule = FALSE)
bad <- vapply(res, function(z) inherits(z, "try-error") || is.null(z), logical(1))
if (any(bad)) stop(sum(bad), " replicates failed: ", as.character(res[bad][[1]]))
saveRDS(list(arms = ARMS, kappa_c = KAPPA_C, nrep = NREP, n = n, p = p,
             rho = rho, jobs = res), out_path, compress = "xz")

RULES <- c("screen, reglage unique", "screen, 9 reglages minimum",
           "quantile SIS tau=.95")
fm <- function(v, k = 3) formatC(v, format = "f", digits = k)
cat(sprintf("\nn=%d p=%d rho=%.2f, %d replications\n\n", n, p, rho, NREP))
for (i in seq_along(ARMS)) {
  w <- res[vapply(res, function(z) z$arm == i, logical(1))]
  sdl <- mean(vapply(w, function(z) z$sd_log_ell, numeric(1)))
  cat(sprintf("=== %s   sd(log ell) = %s\n", ARMS[[i]]$tag,
              if (is.na(sdl)) "-" else fm(sdl)))
  cat(sprintf("%-28s %s %11s %11s %11s\n", "regle",
              paste(sprintf("%8s", paste0("S", DS)), collapse = ""),
              "top4 gamma", "top4 ech.", "top24 ech."))
  for (k in seq_along(RULES)) {
    rm_ <- vapply(w, function(z) z$rmax[k], numeric(1))
    cat(sprintf("%-28s %s %11s %11s %11s\n", RULES[k],
      paste(sprintf("%8s", fm(vapply(DS, function(d) mean(rm_ <= d), numeric(1)))),
            collapse = ""),
      fm(mean(vapply(w, function(z) z$t4g[k], numeric(1))), 2),
      fm(mean(vapply(w, function(z) z$t4s[k], numeric(1))), 2),
      fm(mean(vapply(w, function(z) z$t24s[k], numeric(1))), 2)))
  }
  cat("\n")
}
