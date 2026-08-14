## Follow-up to code/R/test_scale_cluster.R, addressing its two stated caveats.
##
##  (1) Variance matching.  The first run set kappa_c = 2*kappa from the nominal
##      Var(log ell) = kappa^2 s/12, but the AR(1) correlation of the scale block
##      inflates the independent sum beyond that nominal value, so the cluster
##      arm was ~17% less noisy than the s=4 independent baseline it was meant to
##      match.  Here kappa_c is calibrated on the *measured* sd of the baseline:
##      log ell = kappa_c (Phi(F) - 1/2) has sd = kappa_c/sqrt(12) exactly, for
##      every s and lambda, so kappa_c = sd_target * sqrt(12).
##
##  (2) lambda sweep.  lambda fixes how well each observed proxy tracks the
##      latent nuisance, i.e. how attractive each decoy is on its own; it is a
##      second dial and was pinned at 0.9 without justification.  Sweeping it
##      says whether the effect survives imperfect proxies.
##
## usage: KAPPA=0.20 Rscript code/R/test_scale_cluster_sweep.R OUT.rds [NREP] [CORES]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R")
source("code/R/qa_sis.R")

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1L) args[1L] else "results/scale_cluster/sweep.rds"
NREP  <- if (length(args) >= 2L) as.integer(args[2L]) else 40L
CORES <- if (length(args) >= 3L) as.integer(args[3L]) else 90L
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
KAP <- as.numeric(Sys.getenv("KAPPA", "0.20"))
TAU <- 0.95
n <- 2000L; p <- 1000L; rho <- 0.25
SS <- c(4L, 8L, 16L, 32L)
LAMBDAS <- c(0.5, 0.7, 0.9)

draw <- function(seed, s, design, kappa_c, lambda) {
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
  log_ell <- if (design == "cluster") kappa_c * (pnorm(f) - 0.5)
             else KAP * rowSums(u[, scale_idx, drop = FALSE] - 0.5)
  y <- V^(-gamma) * exp(-V / 2) * exp(log_ell)
  stopifnot(all(is.finite(y)), all(y > 0))
  list(z = z, y = y, scale_idx = scale_idx, sd_log_ell = sd(log_ell))
}

## --- calibration: measured sd of the s=4 independent baseline ---------------
cal <- vapply(1:20, function(r)
  draw(31000019L + r * 401L, 4L, "indep", NA_real_, NA_real_)$sd_log_ell,
  numeric(1))
SD_TARGET <- mean(cal)
KAPPA_C <- SD_TARGET * sqrt(12)
cat(sprintf("calibration : sd(log ell) de la reference indep s=4 = %.4f (+- %.4f) -> kappa_c = %.4f\n",
            SD_TARGET, sd(cal), KAPPA_C))

grid <- rbind(
  data.frame(s = SS, design = "indep", lambda = NA_real_,
             kappa_c = NA_real_, stringsAsFactors = FALSE),
  do.call(rbind, lapply(LAMBDAS, function(l)
    data.frame(s = SS, design = "cluster", lambda = l, kappa_c = KAPPA_C,
               stringsAsFactors = FALSE))))
jobs <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
  data.frame(i = i, r = seq_len(NREP))))

one <- function(ix) {
  g <- grid[jobs$i[ix], ]; r <- jobs$r[ix]
  seed <- 31000019L + jobs$i[ix] * 100019L + r * 307L
  d <- draw(seed, g$s, g$design, g$kappa_c, g$lambda)
  h <- n^(-BSTAR) / 2; alpha <- n^(-ASTAR)
  sc <- vapply(seq_len(p), function(j)
    score_coordinate_cpp(d$z[, j], d$y, h, alpha, EPS), numeric(5))["score", ]
  o_t <- order(sc, seq_len(p), na.last = TRUE)
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (n + 1))
  o_q <- order(-qa_sis_scores(uh, d$y, tau = TAU), seq_len(p))
  list(s = g$s, design = g$design, lambda = g$lambda, replicate = r, seed = seed,
       sd_log_ell = d$sd_log_ell, rt = match(1:4, o_t), rq = match(1:4, o_q),
       t4t = o_t[1:4], t4q = o_q[1:4], scale_idx = d$scale_idx)
}

res <- parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = CORES,
                          mc.preschedule = FALSE)
bad <- vapply(res, function(z) inherits(z, "try-error") || is.null(z), logical(1))
if (any(bad)) stop(sum(bad), " replicates failed")
saveRDS(list(grid = grid, nrep = NREP, kappa = KAP, kappa_c = KAPPA_C,
             sd_target = SD_TARGET, n = n, p = p, rho = rho, tau = TAU,
             jobs = res), out_path, compress = "xz")

DS <- c(4, 10, 20, 30, 50, 100)
cat(sprintf("\nM2 : nombre de covariables d'echelle, variances appariees\n"))
cat(sprintf("n=%d p=%d rho=%.2f, %d replications, kappa=%.2f, kappa_c=%.3f, tau=%.2f\n\n",
            n, p, rho, NREP, KAP, KAPPA_C, TAU))
cat(sprintf("%-8s %6s %3s %8s | %-16s %s | %s\n", "design", "lambda", "s",
            "sd(logl)", "methode",
            paste(sprintf("%7s", paste0("S", DS)), collapse = ""),
            "top4 g/ech"))
for (i in seq_len(nrow(grid))) {
  w <- res[vapply(res, function(z) z$s == grid$s[i] &&
                    z$design == grid$design[i] &&
                    identical(z$lambda, grid$lambda[i]), logical(1))]
  if (!length(w)) next
  sdl <- mean(vapply(w, function(z) z$sd_log_ell, numeric(1)))
  first <- TRUE
  for (meth in c("Tail-index SIS", "Quantile SIS .95")) {
    key <- if (meth == "Tail-index SIS") "rt" else "rq"
    t4 <- if (meth == "Tail-index SIS") "t4t" else "t4q"
    rmax <- vapply(w, function(z) max(z[[key]]), numeric(1))
    g4 <- mean(vapply(w, function(z) sum(z[[t4]] %in% 1:4), numeric(1)))
    s4 <- mean(vapply(w, function(z) sum(z[[t4]] %in% z$scale_idx), numeric(1)))
    cat(sprintf("%-8s %6s %3s %8s | %-16s %s | %.2f / %.2f\n",
                if (first) grid$design[i] else "",
                if (first && !is.na(grid$lambda[i])) sprintf("%.1f", grid$lambda[i]) else "",
                if (first) grid$s[i] else "",
                if (first) sprintf("%.3f", sdl) else "",
                meth,
                paste(sprintf("%7.2f", vapply(DS, function(d) mean(rmax <= d),
                                              numeric(1))), collapse = ""),
                g4, s4))
    first <- FALSE
  }
}
cat("\nWROTE", out_path, "\n")
