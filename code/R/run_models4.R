## Validation de la suite de modèles figée (code/R/generate4.R).
## Vérifie que la nuisance de M2 a bien la magnitude visée et indépendante de
## rho, puis mesure les trois règles sur les quatre modèles.
## usage: GRID_N=2000 GRID_P=1000 GRID_RHO=0.25 \
##        Rscript code/R/run_models4.R OUT.rds [NREP] [CORES]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate4.R"); source("code/R/qa_sis.R")

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1L) args[1L] else "results/models4/raw.rds"
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
D1 <- 25L
MODELS <- c("M1", "M2", "M3", "M4")
RULES <- c("screen seul", "screen 9 minimum", "quantile .95",
           "pipeline + 9 min")

score_vec <- function(z, y, a, b) {
  h <- n^(-b) / 2; alpha <- n^(-a)
  vapply(seq_len(p), function(j)
    score_coordinate_cpp(z[, j], y, h, alpha, EPS), numeric(5))["score", ]
}

jobs <- do.call(rbind, lapply(seq_along(MODELS), function(i)
  data.frame(arm = i, r = seq_len(NREP))))

one <- function(ix) {
  m <- MODELS[jobs$arm[ix]]; r <- jobs$r[ix]
  seed <- 111000019L + jobs$arm[ix] * 100003L + r * 307L
  d <- simulate_dataset4(n, p, rho, m, seed)
  base <- score_vec(d$z, d$y, ASTAR, BSTAR)
  R <- matrix(NA_real_, p, length(AGRID) * length(BGRID)); ii <- 0L
  for (a in AGRID) for (b in BGRID) {
    ii <- ii + 1L
    R[, ii] <- rank(score_vec(d$z, d$y, a, b), ties.method = "first",
                    na.last = TRUE)
  }
  amin <- apply(R, 1L, min)
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (n + 1))
  ord_q <- order(-qa_sis_scores(uh, d$y, tau = 0.95), seq_len(p))
  surv <- ord_q[seq_len(D1)]
  r_pipe <- if (!all(A_GAMMA %in% surv)) {
    Inf
  } else {
    max(match(A_GAMMA, surv[order(amin[surv], surv)]))
  }
  ords <- list(order(base, seq_len(p), na.last = TRUE),
               order(amin, seq_len(p)), ord_q)
  list(model = m, replicate = r, seed = seed,
       gmin = min(d$gamma), gmax = max(d$gamma),
       sd_log_ell = if (m == "M2")
         sd(M2_KAPPA_C * (pnorm(d$f) - 0.5)) else NA_real_,
       rmax = c(vapply(ords, function(o) max(match(A_GAMMA, o)), numeric(1)),
                r_pipe),
       t4s = c(vapply(ords, function(o) sum(o[1:4] %in% A_SCALE), numeric(1)),
               NA_real_),
       t24s = c(vapply(ords, function(o) sum(o[1:24] %in% A_SCALE), numeric(1)),
                NA_real_))
}

res <- parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = CORES,
                          mc.preschedule = FALSE)
bad <- vapply(res, function(z) inherits(z, "try-error") || is.null(z), logical(1))
if (any(bad)) stop(sum(bad), " replicates failed: ", as.character(res[bad][[1]]))
saveRDS(list(models = MODELS, n = n, p = p, rho = rho, nrep = NREP,
             lambda = M2_LAMBDA, s = M2_S, kappa_c = M2_KAPPA_C, jobs = res),
        out_path, compress = "xz")

fm <- function(v, k = 3) formatC(v, format = "f", digits = k)
cat(sprintf("\nSuite figee : lambda=%.2f, s=%d, kappa_c=%.4f -> sd(log ell) visee %.4f\n",
            M2_LAMBDA, M2_S, M2_KAPPA_C, M2_KAPPA_C / sqrt(12)))
cat(sprintf("n=%d p=%d rho=%.2f, %d replications, A_scale = {%d..%d}\n\n",
            n, p, rho, NREP, min(A_SCALE), max(A_SCALE)))
for (i in seq_along(MODELS)) {
  w <- res[vapply(res, function(z) z$model == MODELS[i], logical(1))]
  sdl <- mean(vapply(w, function(z) z$sd_log_ell, numeric(1)))
  cat(sprintf("=== %s   gamma dans (%.3f, %.3f)%s\n", MODELS[i],
              mean(vapply(w, function(z) z$gmin, numeric(1))),
              mean(vapply(w, function(z) z$gmax, numeric(1))),
              if (is.na(sdl)) "" else sprintf(", sd(log ell) = %.4f", sdl)))
  cat(sprintf("%-22s %8s %8s %8s %11s %11s\n", "regle", "Sure-4", "Sure-10",
              "Sure-20", "top4 ech.", "top24 ech."))
  for (k in seq_along(RULES)) {
    v <- vapply(w, function(z) z$rmax[k], numeric(1))
    cat(sprintf("%-22s %8s %8s %8s %11s %11s\n", RULES[k],
                fm(mean(v <= 4)), fm(mean(v <= 10)), fm(mean(v <= 20)),
                if (k == 4) "-" else fm(mean(vapply(w, function(z) z$t4s[k], numeric(1))), 2),
                if (k == 4) "-" else fm(mean(vapply(w, function(z) z$t24s[k], numeric(1))), 2)))
  }
  cat("\n")
}
