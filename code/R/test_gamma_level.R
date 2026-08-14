## Raising the leading constant of gamma in M3 and M4 to 1.
##
##   M3 publie : gamma = 0.50 * exp(-0.50 * s4)      -> gamma in (0.068, 0.50)
##   M3 modifie: gamma = 1.00 * exp(-0.50 * s4)      -> gamma in (0.135, 1.00)
##   M4 publie : gamma = 0.55 * exp(-0.35 s4 - 0.20 u1u2 - 0.15 u3u4)
##   M4 modifie: gamma = 1.00 * exp(-0.35 s4 - 0.20 u1u2 - 0.15 u3u4)
##
## The RELATIVE contrast along a fibre is unchanged (a constant multiplies the
## envelope of active and inactive coordinates alike, so Delta and the Hill
## sampling noise both scale with gamma).  What changes is the absolute tail
## heaviness: log Y = gamma log(1/V) + log ell, and log ell is bounded with a
## fixed variance, so doubling gamma halves the nuisance relative to the signal.
##
## usage: KAPPA=0.20 GRID_N=2000 GRID_P=1000 GRID_RHO=0 \
##        Rscript code/R/test_gamma_level.R OUT.rds [NREP] [CORES]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R"); source("code/R/qa_sis.R")

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1L) args[1L] else "results/gamma_level/raw.rds"
NREP  <- if (length(args) >= 2L) as.integer(args[2L]) else 40L
CORES <- if (length(args) >= 3L) as.integer(args[3L]) else 92L
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.05, 0.10, 0.15)
n <- as.integer(Sys.getenv("GRID_N", "2000"))
p <- as.integer(Sys.getenv("GRID_P", "1000"))
rho <- as.numeric(Sys.getenv("GRID_RHO", "0"))
D1 <- 25L
DS <- c(4, 10, 20)
RULES <- c("screen seul", "screen 9 minimum", "quantile .95",
           "pipeline + 9 min")

gam <- function(u, model, lead) {
  s4 <- rowSums(u[, 1:4, drop = FALSE])
  if (model == "M3") lead * exp(-0.5 * s4)
  else lead * exp(-0.35 * s4 - 0.20 * u[, 1] * u[, 2] - 0.15 * u[, 3] * u[, 4])
}

ARMS <- list(
  list(model = "M3", lead = 0.50, tag = "M3 publie   (0.50)"),
  list(model = "M3", lead = 1.00, tag = "M3 modifie  (1.00)"),
  list(model = "M4", lead = 0.55, tag = "M4 publie   (0.55)"),
  list(model = "M4", lead = 1.00, tag = "M4 modifie  (1.00)"))

score_vec <- function(z, y, a, b) {
  h <- n^(-b) / 2; alpha <- n^(-a)
  vapply(seq_len(p), function(j)
    score_coordinate_cpp(z[, j], y, h, alpha, EPS), numeric(5))["score", ]
}

jobs <- do.call(rbind, lapply(seq_along(ARMS), function(i)
  data.frame(arm = i, r = seq_len(NREP))))

one <- function(ix) {
  A <- ARMS[[jobs$arm[ix]]]; r <- jobs$r[ix]
  seed <- 91000003L + jobs$arm[ix] * 100003L + r * 307L
  set.seed(seed)
  z <- ar1_gaussian(n, p, rho)
  u <- population_uniform_matrix(z)
  gamma <- gam(u, A$model, A$lead)
  V <- runif(n)
  ## published slowly varying factor of M1/M3/M4, unchanged
  ell <- 1 / (1 + exp(sv_l1(u) - pmin(1 / V, 700)))
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
       gamma_min = min(gamma), gamma_max = max(gamma),
       rmax = c(max(match(1:4, order(base, seq_len(p), na.last = TRUE))),
                max(match(1:4, order(amin, seq_len(p)))),
                max(match(1:4, ord_q)),
                r_pipe))
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
  w <- res[vapply(res, function(z) z$arm == i, logical(1))]
  cat(sprintf("=== %s   gamma dans (%.3f, %.3f)\n", ARMS[[i]]$tag,
              mean(vapply(w, function(z) z$gamma_min, numeric(1))),
              mean(vapply(w, function(z) z$gamma_max, numeric(1)))))
  cat(sprintf("%-22s %8s %8s %8s\n", "regle", "Sure-4", "Sure-10", "Sure-20"))
  for (k in seq_along(RULES)) {
    v <- vapply(w, function(z) z$rmax[k], numeric(1))
    cat(sprintf("%-22s %8s %8s %8s\n", RULES[k],
                fm(mean(v <= 4)), fm(mean(v <= 10)), fm(mean(v <= 20))))
  }
  cat("\n")
}
