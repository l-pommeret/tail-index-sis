## M2 with 20 scale-active covariates instead of 4.
##
## A_gamma = {1,2,3,4} drives the tail index; A_scale = {5,...,4+s} moves only
## the finite conditional quantiles, through log ell = kappa * sum(u_j - 1/2).
## Var(log ell) = kappa^2 s / 12, so raising s at fixed kappa also raises the
## nuisance magnitude; two arms separate the two effects:
##   naif    s = 20, kappa = 0.20 unchanged  -> nuisance sd multiplied by sqrt(5)
##   apparie s = 20, kappa = 0.20*sqrt(4/20) -> nuisance sd held at its s=4 value
## plus the s = 4 reference, which is the published M2.
##
## Screens compared: the proposed score at its single published tuning, the same
## score aggregated by MINIMUM RANK over the 3x3 tuning grid, and quantile SIS
## at tau = 0.95.
##
## usage: KAPPA=0.20 Rscript code/R/test_m2_scale20.R OUT.rds [NREP] [CORES]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R"); source("code/R/qa_sis.R")

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1L) args[1L] else "results/m2_scale20/raw.rds"
NREP  <- if (length(args) >= 2L) as.integer(args[2L]) else 40L
CORES <- if (length(args) >= 3L) as.integer(args[3L]) else 92L
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
KAP <- as.numeric(Sys.getenv("KAPPA", "0.20"))
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.05, 0.10, 0.15)   # the 9 tunings
n <- 2000L; p <- 1000L; rho <- 0.25
DS <- c(4, 10, 20, 30, 50)

ARMS <- list(
  list(tag = "s=4 (M2 publie)",   s = 4L,  kappa = KAP),
  list(tag = "s=20 kappa inchange", s = 20L, kappa = KAP),
  list(tag = "s=20 variance appariee", s = 20L, kappa = KAP * sqrt(4 / 20)))

draw <- function(seed, s, kappa) {
  set.seed(seed)
  z <- ar1_gaussian(n, p, rho)
  u <- population_uniform_matrix(z)
  scale_idx <- 5:(4 + s)
  gamma <- 0.5 * exp(-rowSums(u[, 1:4, drop = FALSE]))
  V <- runif(n)
  log_ell <- kappa * rowSums(u[, scale_idx, drop = FALSE] - 0.5)
  y <- V^(-gamma) * exp(-V / 2) * exp(log_ell)
  stopifnot(all(is.finite(y)), all(y > 0))
  list(z = z, y = y, scale_idx = scale_idx, sd_log_ell = sd(log_ell))
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
  seed <- 41000003L + jobs$arm[ix] * 100003L + r * 307L
  d <- draw(seed, A$s, A$kappa)
  base <- score_vec(d$z, d$y, ASTAR, BSTAR)
  R <- matrix(NA_real_, p, length(AGRID) * length(BGRID)); ii <- 0L
  for (a in AGRID) for (b in BGRID) {
    ii <- ii + 1L
    R[, ii] <- rank(score_vec(d$z, d$y, a, b), ties.method = "first",
                    na.last = TRUE)
  }
  agg <- apply(R, 1L, min)                       # minimum rank over 9 tunings
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (n + 1))
  ords <- list(order(base, seq_len(p), na.last = TRUE),
               order(agg, seq_len(p)),
               order(-qa_sis_scores(uh, d$y, tau = 0.95), seq_len(p)))
  list(arm = jobs$arm[ix], replicate = r, seed = seed, s = A$s,
       sd_log_ell = d$sd_log_ell,
       rmax = vapply(ords, function(o) max(match(1:4, o)), numeric(1)),
       t4g = vapply(ords, function(o) sum(o[1:4] %in% 1:4), numeric(1)),
       t4s = vapply(ords, function(o) sum(o[1:4] %in% d$scale_idx), numeric(1)),
       top24s = vapply(ords, function(o) sum(o[1:24] %in% d$scale_idx), numeric(1)))
}

res <- parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = CORES,
                          mc.preschedule = FALSE)
bad <- vapply(res, function(z) inherits(z, "try-error") || is.null(z), logical(1))
if (any(bad)) stop(sum(bad), " replicates failed: ", as.character(res[bad][[1]]))
saveRDS(list(arms = ARMS, nrep = NREP, n = n, p = p, rho = rho, jobs = res),
        out_path, compress = "xz")

RULES <- c("screen, reglage unique", "screen, 9 reglages minimum",
           "quantile SIS tau=.95")
fm <- function(v, k = 3) formatC(v, format = "f", digits = k)
cat(sprintf("\nM2 avec s covariables d'echelle : n=%d p=%d rho=%.2f, %d replications\n\n",
            n, p, rho, NREP))
for (i in seq_along(ARMS)) {
  w <- res[vapply(res, function(z) z$arm == i, logical(1))]
  cat(sprintf("=== %s : A_scale = {5..%d}, kappa = %.4f, sd(log ell) = %.3f\n",
              ARMS[[i]]$tag, 4 + ARMS[[i]]$s, ARMS[[i]]$kappa,
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
      fm(mean(vapply(w, function(z) z$top24s[k], numeric(1))), 2)))
  }
  cat("\n")
}
