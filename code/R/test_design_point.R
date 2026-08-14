## Evaluate the screens at an arbitrary design point, generating fresh data
## rather than replaying campaign cells (needed for n outside {1000,2000,5000}).
## Rules: proposed score at its published tuning; the same aggregated by median
## and by minimum rank over the 3x3 tuning grid; quantile SIS at tau = 0.95;
## and the two-stage pipeline behind the minimum-rank aggregation.
##
## Seed stream 71000003 + arm*100003 + r*307, disjoint from every other stream.
##
## usage: KAPPA=0.20 GRID_N=4000 GRID_P=1000 GRID_RHO=0 \
##        Rscript code/R/test_design_point.R OUT.rds MODELS [NREP] [CORES]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R"); source("code/R/qa_sis.R")

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1L) args[1L] else "results/design_point/raw.rds"
MODELS <- if (length(args) >= 2L) {
  strsplit(args[2L], ",")[[1]]
} else {
  c("M1", "M2", "M3", "M4")
}
NREP  <- if (length(args) >= 3L) as.integer(args[3L]) else 40L
CORES <- if (length(args) >= 4L) as.integer(args[4L]) else 92L
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.05, 0.10, 0.15)
n <- as.integer(Sys.getenv("GRID_N", "2000"))
p <- as.integer(Sys.getenv("GRID_P", "1000"))
rho <- as.numeric(Sys.getenv("GRID_RHO", "0.25"))
D1 <- 25L
DS <- c(4, 10, 20, 30, 50)
RULES <- c("screen seul", "screen 9 regl. mediane", "screen 9 regl. minimum",
           "quantile SIS .95", "pipeline + 9 regl. min")

score_vec <- function(z, y, a, b) {
  h <- n^(-b) / 2; alpha <- n^(-a)
  vapply(seq_len(p), function(j)
    score_coordinate_cpp(z[, j], y, h, alpha, EPS), numeric(5))["score", ]
}

jobs <- do.call(rbind, lapply(seq_along(MODELS), function(i)
  data.frame(arm = i, r = seq_len(NREP))))

one <- function(ix) {
  m <- MODELS[jobs$arm[ix]]; r <- jobs$r[ix]
  seed <- 71000003L + jobs$arm[ix] * 100003L + r * 307L
  d <- simulate_dataset3(n, p, rho, m, seed)
  base <- score_vec(d$z, d$y, ASTAR, BSTAR)
  R <- matrix(NA_real_, p, length(AGRID) * length(BGRID)); ii <- 0L
  for (a in AGRID) for (b in BGRID) {
    ii <- ii + 1L
    R[, ii] <- rank(score_vec(d$z, d$y, a, b), ties.method = "first",
                    na.last = TRUE)
  }
  amed <- apply(R, 1L, median); amin <- apply(R, 1L, min)
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (n + 1))
  ord_q <- order(-qa_sis_scores(uh, d$y, tau = 0.95), seq_len(p))
  surv <- ord_q[seq_len(D1)]
  r_pipe <- if (!all(1:4 %in% surv)) {
    Inf
  } else {
    max(match(1:4, surv[order(amin[surv], surv)]))
  }
  rmax <- c(max(match(1:4, order(base, seq_len(p), na.last = TRUE))),
            max(match(1:4, order(amed, seq_len(p)))),
            max(match(1:4, order(amin, seq_len(p)))),
            max(match(1:4, ord_q)),
            r_pipe)
  list(model = m, replicate = r, seed = seed, rmax = rmax)
}

res <- parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = CORES,
                          mc.preschedule = FALSE)
bad <- vapply(res, function(z) inherits(z, "try-error") || is.null(z), logical(1))
if (any(bad)) stop(sum(bad), " replicates failed: ", as.character(res[bad][[1]]))
saveRDS(list(models = MODELS, n = n, p = p, rho = rho, nrep = NREP, jobs = res),
        out_path, compress = "xz")

fm <- function(v) formatC(v, format = "f", digits = 3)
cat(sprintf("\nn=%d p=%d rho=%.2f, %d replications\n\n", n, p, rho, NREP))
cat(sprintf("%-24s %s\n", "regle",
            paste(sprintf("%-17s", MODELS), collapse = "")))
for (k in seq_along(RULES)) {
  vals <- vapply(MODELS, function(m) {
    w <- res[vapply(res, function(z) z$model == m, logical(1))]
    v <- vapply(w, function(z) z$rmax[k], numeric(1))
    sprintf("%s/%s", fm(mean(v <= 4)), fm(mean(v <= 20)))
  }, character(1))
  cat(sprintf("%-24s %s\n", RULES[k],
              paste(sprintf("%-17s", vals), collapse = "")))
}
cat("\n(Sure-4 / Sure-20)\n")
