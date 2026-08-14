## Raising M4's decay rate so that its separation matches M1 and M2.
##
## For gamma = c exp(-beta*s4 - interactions), conditioning on u_j = t and
## maximising over the other three active coordinates puts them at 0, where both
## interaction terms vanish.  The upper envelope is therefore exactly
## c exp(-beta t), so M4's effective decay rate is 0.35 against 1.00 for M1/M2,
## and its population separation is
##   Delta/c = 1 - {e^(-beta eps) - e^(-beta(1-eps))} / {beta (1-2 eps)}
##           = 0.157 for beta=0.35 against 0.373 for beta=1.
## Raising beta to 1 should bring M4 to M1's level while keeping its
## interactions, which live in the interior of the surface and not on the
## envelope.
##
## Arms:
##   M1                reference, beta = 1
##   M4 publie         beta = 0.35, interactions 0.20 / 0.15
##   M4 beta=1         beta raised, interactions left as published
##   M4 exposant x2.86 whole exponent scaled by 1/0.35, so the interactions keep
##                     their relative weight (same envelope, same separation)
##
## usage: KAPPA=0.20 GRID_N=2000 GRID_P=1000 GRID_RHO=0.25 \
##        Rscript code/R/test_m4_beta.R OUT.rds [NREP] [CORES]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R"); source("code/R/qa_sis.R")

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1L) args[1L] else "results/m4_beta/raw.rds"
NREP  <- if (length(args) >= 2L) as.integer(args[2L]) else 40L
CORES <- if (length(args) >= 3L) as.integer(args[3L]) else 92L
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.05, 0.10, 0.15)
n <- as.integer(Sys.getenv("GRID_N", "2000"))
p <- as.integer(Sys.getenv("GRID_P", "1000"))
rho <- as.numeric(Sys.getenv("GRID_RHO", "0.25"))
D1 <- 25L; DS <- c(4, 10, 20)
RULES <- c("screen seul", "screen 9 minimum", "quantile .95", "pipeline + 9 min")

sep <- function(b) 1 - (exp(-b * EPS) - exp(-b * (1 - EPS))) / (b * (1 - 2 * EPS))

## beta grid overridable: BETAS="0.35,0.50,0.60,0.70,0.85,1.00"
BETAS <- as.numeric(strsplit(Sys.getenv("BETAS", "0.35,1.00"), ",")[[1]])
ARMS <- c(
  list(list(tag = "M1 (reference, beta=1)", kind = "M1", beta = 1.00,
            i12 = 0, i34 = 0, c = 0.50)),
  lapply(BETAS, function(b)
    list(tag = sprintf("M4 beta=%.2f (interactions publiees)", b), kind = "M4",
         beta = b, i12 = 0.20, i34 = 0.15, c = 0.55)))

gam <- function(u, A) {
  s4 <- rowSums(u[, 1:4, drop = FALSE])
  if (A$kind == "M1") A$c * exp(-A$beta * s4)
  else A$c * exp(-A$beta * s4 - A$i12 * u[, 1] * u[, 2] - A$i34 * u[, 3] * u[, 4])
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
  seed <- 101000003L + jobs$arm[ix] * 100003L + r * 307L
  set.seed(seed)
  z <- ar1_gaussian(n, p, rho)
  u <- population_uniform_matrix(z)
  gamma <- gam(u, A)
  V <- runif(n)
  ell <- 1 / (1 + exp(sv_l1(u) - pmin(1 / V, 700)))   # published l1, unchanged
  y <- V^(-gamma) * ell
  stopifnot(all(is.finite(y)), all(y > 0))
  base <- score_vec(z, y, ASTAR, BSTAR)
  R <- matrix(NA_real_, p, length(AGRID) * length(BGRID)); ii <- 0L
  for (a in AGRID) for (b in BGRID) {
    ii <- ii + 1L
    R[, ii] <- rank(score_vec(z, y, a, b), ties.method = "first", na.last = TRUE)
  }
  amin <- apply(R, 1L, min)
  uh <- apply(z, 2L, function(x) rank(x, ties.method = "average") / (n + 1))
  ord_q <- order(-qa_sis_scores(uh, y, tau = 0.95), seq_len(p))
  surv <- ord_q[seq_len(D1)]
  r_pipe <- if (!all(1:4 %in% surv)) {
    Inf
  } else {
    max(match(1:4, surv[order(amin[surv], surv)]))
  }
  list(arm = jobs$arm[ix], replicate = r, seed = seed,
       gmin = min(gamma), gmax = max(gamma),
       rmax = c(max(match(1:4, order(base, seq_len(p), na.last = TRUE))),
                max(match(1:4, order(amin, seq_len(p)))),
                max(match(1:4, ord_q)), r_pipe))
}

res <- parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = CORES,
                          mc.preschedule = FALSE)
bad <- vapply(res, function(z) inherits(z, "try-error") || is.null(z), logical(1))
if (any(bad)) stop(sum(bad), " replicates failed: ", as.character(res[bad][[1]]))
saveRDS(list(arms = ARMS, n = n, p = p, rho = rho, nrep = NREP, jobs = res),
        out_path, compress = "xz")

fm <- function(v) formatC(v, format = "f", digits = 3)
cat(sprintf("\nn=%d p=%d rho=%.2f, %d replications\n\n", n, p, rho, NREP))
for (i in seq_along(ARMS)) {
  A <- ARMS[[i]]
  w <- res[vapply(res, function(z) z$arm == i, logical(1))]
  cat(sprintf("=== %s\n    beta = %.2f, separation population Delta/c = %.3f, gamma dans (%.3f, %.3f)\n",
              A$tag, A$beta, sep(A$beta),
              mean(vapply(w, function(z) z$gmin, numeric(1))),
              mean(vapply(w, function(z) z$gmax, numeric(1)))))
  cat(sprintf("%-22s %8s %8s %8s\n", "regle", "Sure-4", "Sure-10", "Sure-20"))
  for (k in seq_along(RULES)) {
    v <- vapply(w, function(z) z$rmax[k], numeric(1))
    cat(sprintf("%-22s %8s %8s %8s\n", RULES[k],
                fm(mean(v <= 4)), fm(mean(v <= 10)), fm(mean(v <= 20))))
  }
  cat("\n")
}
