## M2 with 20 scale covariates carried by a correlated cluster.
##
## The independent-sum version of M2 ties two things to the same coefficient:
## how attractive each scale covariate is to a quantile screen (proportional to
## kappa) and how much the nuisance damages the local Hill statistic
## (proportional to kappa^2 s).  Raising s at fixed kappa therefore kills both
## screens, and matching the variance by kappa ~ 1/sqrt(s) makes each decoy too
## weak to fool anything -- measured in code/R/test_m2_scale20.R.
##
## The cluster version breaks the tie: one latent factor F carries the whole
## nuisance, and the s observed scale covariates are proxies correlated with it,
##   z_j = lambda F + sqrt(1-lambda^2) eps_j,   j in A_scale
##   log ell = kappa_c (Phi(F) - 1/2)
## so Var(log ell) = kappa_c^2/12 whatever s and lambda, while each proxy keeps
## marginal correlation lambda with the nuisance.  kappa_c is calibrated on the
## MEASURED sd of the published s=4 arm, so all arms are equally noisy.
##
## Screens: proposed score at its published tuning, the same aggregated by
## MINIMUM RANK over the 3x3 tuning grid, and quantile SIS at tau = 0.95.
##
## usage: KAPPA=0.20 Rscript code/R/test_m2_cluster20.R OUT.rds [NREP] [CORES]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R"); source("code/R/qa_sis.R")

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1L) args[1L] else "results/m2_cluster20/raw.rds"
NREP  <- if (length(args) >= 2L) as.integer(args[2L]) else 40L
CORES <- if (length(args) >= 3L) as.integer(args[3L]) else 92L
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
KAP <- as.numeric(Sys.getenv("KAPPA", "0.20"))
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.05, 0.10, 0.15)
n <- 2000L; p <- 1000L; rho <- 0.25
DS <- c(4, 10, 20, 30, 50)
LAMBDAS <- c(0.5, 0.7, 0.9)
S_CLUSTER <- 20L

draw <- function(seed, s, design, kappa, lambda = NA_real_) {
  set.seed(seed)
  z <- ar1_gaussian(n, p, rho)
  scale_idx <- 5:(4 + s)
  f <- NULL
  if (design == "cluster") {
    f <- rnorm(n)
    for (j in scale_idx)
      z[, j] <- lambda * f + sqrt(1 - lambda^2) * rnorm(n)
  }
  u <- population_uniform_matrix(z)
  gamma <- 0.5 * exp(-rowSums(u[, 1:4, drop = FALSE]))
  V <- runif(n)
  log_ell <- if (design == "cluster") kappa * (pnorm(f) - 0.5)
             else kappa * rowSums(u[, scale_idx, drop = FALSE] - 0.5)
  y <- V^(-gamma) * exp(-V / 2) * exp(log_ell)
  stopifnot(all(is.finite(y)), all(y > 0))
  list(z = z, y = y, scale_idx = scale_idx, sd_log_ell = sd(log_ell))
}

## calibrate kappa_c on the measured sd of the published s=4 arm
cal <- vapply(1:20, function(r)
  draw(51000017L + r * 401L, 4L, "indep", KAP)$sd_log_ell, numeric(1))
SD_TARGET <- mean(cal); KAPPA_C <- SD_TARGET * sqrt(12)
cat(sprintf("calibration : sd(log ell) de M2 publie = %.4f -> kappa_c = %.4f\n",
            SD_TARGET, KAPPA_C))

ARMS <- c(list(list(tag = "s=4 independant (M2 publie)", s = 4L,
                    design = "indep", kappa = KAP, lambda = NA_real_)),
          lapply(LAMBDAS, function(l)
            list(tag = sprintf("s=20 amas correle, lambda=%.1f", l),
                 s = S_CLUSTER, design = "cluster", kappa = KAPPA_C, lambda = l)))

score_vec <- function(z, y, a, b) {
  h <- n^(-b) / 2; alpha <- n^(-a)
  vapply(seq_len(p), function(j)
    score_coordinate_cpp(z[, j], y, h, alpha, EPS), numeric(5))["score", ]
}

jobs <- do.call(rbind, lapply(seq_along(ARMS), function(i)
  data.frame(arm = i, r = seq_len(NREP))))

one <- function(ix) {
  A <- ARMS[[jobs$arm[ix]]]; r <- jobs$r[ix]
  seed <- 51000017L + jobs$arm[ix] * 100003L + r * 307L
  d <- draw(seed, A$s, A$design, A$kappa, A$lambda)
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
saveRDS(list(arms = ARMS, kappa_c = KAPPA_C, sd_target = SD_TARGET,
             nrep = NREP, n = n, p = p, rho = rho, jobs = res),
        out_path, compress = "xz")

RULES <- c("screen, reglage unique", "screen, 9 reglages minimum",
           "quantile SIS tau=.95")
fm <- function(v, k = 3) formatC(v, format = "f", digits = k)
cat(sprintf("\nM2 en amas correle : n=%d p=%d rho=%.2f, %d replications\n\n",
            n, p, rho, NREP))
for (i in seq_along(ARMS)) {
  w <- res[vapply(res, function(z) z$arm == i, logical(1))]
  cat(sprintf("=== %s   sd(log ell) = %.3f\n", ARMS[[i]]$tag,
              mean(vapply(w, function(z) z$sd_log_ell, numeric(1)))))
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
