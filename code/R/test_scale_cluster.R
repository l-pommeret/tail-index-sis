## Can the number of scale-only covariates in M2 be turned into a dial that is
## independent of how much the nuisance hurts the tail-index screen?
##
## M2 carries its scale nuisance as log ell = kappa * sum_{j in A_scale}(u_j-1/2),
## so Var(log ell) = kappa^2 s / 12: adding scale covariates at fixed kappa is
## algebraically the same move as raising kappa, which the calibration pilot
## shows kills both screens at once.  The alternative tested here puts the whole
## nuisance on a single latent factor F, observed through s correlated proxies:
##
##   indep    z_{A_scale} AR(1) as usual,  log ell = kappa * sum (u_j - 1/2)
##   cluster  z_j = lambda F + sqrt(1-lambda^2) eps_j  for j in A_scale,
##            log ell = kappa_c * (Phi(F) - 1/2),  kappa_c = 2 kappa
##
## kappa_c = 2 kappa matches Var(log ell) to the s = 4 independent case for every
## s, so the damage done to the local Hill statistic is held fixed while s grows,
## and each proxy keeps a marginal correlation lambda with the nuisance, so each
## stays individually attractive to a quantile screen.
##
## In the cluster design the scale columns are overwritten after the AR(1) draw,
## so the block's correlation with its two neighbouring columns is by
## construction not rho; everything outside the block is untouched.
##
## usage: KAPPA=0.20 Rscript code/R/test_scale_cluster.R OUT.rds [NREP] [CORES]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R")
source("code/R/qa_sis.R")

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1L) args[1L] else "results/scale_cluster/raw.rds"
NREP  <- if (length(args) >= 2L) as.integer(args[2L]) else 40L
CORES <- if (length(args) >= 3L) as.integer(args[3L]) else 12L
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
KAP <- as.numeric(Sys.getenv("KAPPA", "0.20"))
LAMBDA <- 0.9
TAU <- 0.95
n <- 2000L; p <- 1000L; rho <- 0.25          # the paper's baseline design
SS <- c(4L, 8L, 16L, 32L)
DESIGNS <- c("indep", "cluster")

simulate_variant <- function(n, p, rho, seed, s, design) {
  set.seed(seed)
  z <- ar1_gaussian(n, p, rho)
  scale_idx <- 5:(4 + s)
  if (design == "cluster") {
    f <- rnorm(n)
    for (j in scale_idx)
      z[, j] <- LAMBDA * f + sqrt(1 - LAMBDA^2) * rnorm(n)
  }
  u <- population_uniform_matrix(z)
  gamma <- 0.5 * exp(-rowSums(u[, 1:4, drop = FALSE]))
  V <- runif(n)
  log_ell <- if (design == "cluster")
    2 * KAP * (pnorm(f) - 0.5)
  else
    KAP * rowSums(u[, scale_idx, drop = FALSE] - 0.5)
  y <- V^(-gamma) * exp(-V / 2) * exp(log_ell)
  stopifnot(all(is.finite(y)), all(y > 0))
  list(z = z, y = y, scale_idx = scale_idx, sd_log_ell = sd(log_ell))
}

grid <- expand.grid(s = SS, design = DESIGNS, stringsAsFactors = FALSE)
jobs <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
  data.frame(i = i, r = seq_len(NREP))))

one <- function(ix) {
  g <- grid[jobs$i[ix], ]; r <- jobs$r[ix]
  seed <- 21000011L + jobs$i[ix] * 100003L + r * 307L
  d <- simulate_variant(n, p, rho, seed, g$s, g$design)
  h <- n^(-BSTAR) / 2; alpha <- n^(-ASTAR)
  sc <- vapply(seq_len(p), function(j)
    score_coordinate_cpp(d$z[, j], d$y, h, alpha, EPS), numeric(5))["score", ]
  o_t <- order(sc, seq_len(p), na.last = TRUE)
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (n + 1))
  o_q <- order(-qa_sis_scores(uh, d$y, tau = TAU), seq_len(p))
  list(s = g$s, design = g$design, replicate = r, seed = seed,
       sd_log_ell = d$sd_log_ell,
       rt = match(1:4, o_t), rq = match(1:4, o_q),
       t4t = o_t[1:4], t4q = o_q[1:4], scale_idx = d$scale_idx)
}

res <- parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = CORES,
                          mc.preschedule = FALSE)
bad <- vapply(res, function(z) inherits(z, "try-error") || is.null(z), logical(1))
if (any(bad)) stop(sum(bad), " replicates failed")
saveRDS(list(grid = grid, nrep = NREP, kappa = KAP, lambda = LAMBDA,
             n = n, p = p, rho = rho, tau = TAU, jobs = res),
        out_path, compress = "xz")

DS <- c(4, 10, 20, 30, 50, 100)
cat(sprintf("\nM2 : nombre de covariables d'echelle comme cadran\n"))
cat(sprintf("n=%d p=%d rho=%.2f, %d replications, kappa=%.2f, lambda=%.2f, quantile SIS a tau=%.2f\n\n",
            n, p, rho, NREP, KAP, LAMBDA, TAU))
cat(sprintf("%-8s %3s %8s | %-14s %s | %s\n", "design", "s", "sd(logl)",
            "methode", paste(sprintf("%7s", paste0("S", DS)), collapse = ""),
            "top4: gamma / echelle"))
for (i in seq_len(nrow(grid))) {
  w <- res[vapply(res, function(z) z$s == grid$s[i] &&
                    z$design == grid$design[i], logical(1))]
  sdl <- mean(vapply(w, function(z) z$sd_log_ell, numeric(1)))
  for (meth in c("Tail-index SIS", "Quantile SIS .95")) {
    key <- if (meth == "Tail-index SIS") "rt" else "rq"
    t4 <- if (meth == "Tail-index SIS") "t4t" else "t4q"
    rmax <- vapply(w, function(z) max(z[[key]]), numeric(1))
    g4 <- mean(vapply(w, function(z) sum(z[[t4]] %in% 1:4), numeric(1)))
    s4 <- mean(vapply(w, function(z) sum(z[[t4]] %in% z$scale_idx), numeric(1)))
    cat(sprintf("%-8s %3d %8.3f | %-14s %s | %.2f / %.2f\n",
                if (meth == "Tail-index SIS") grid$design[i] else "",
                if (meth == "Tail-index SIS") grid$s[i] else NA_integer_,
                if (meth == "Tail-index SIS") sdl else NA_real_,
                meth,
                paste(sprintf("%7.2f", vapply(DS, function(d) mean(rmax <= d),
                                              numeric(1))), collapse = ""),
                g4, s4))
  }
}
cat("\nWROTE", out_path, "\n")
