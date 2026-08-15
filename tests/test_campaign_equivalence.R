## Verify that the Draft-5 models A1 and A2 reproduce, at the campaign5
## seeds, exactly the M1 and M3 datasets and rule outputs stored in
## results/campaign5 (so those cells can be relabelled instead of rerun).
## Checks 3 replications in each (model, p) cell against the stored rmax
## of the six campaign5 rules that Draft 5 retains, accounting for the
## campaign5 rule order (screen seul, screen 9 min, YU, q90, q95, q975,
## q99) versus campaign6 (q975 dropped).
suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate5.R"); source("code/R/qa_sis.R")
source("code/R/yoshida_umezu.R")

N <- 2000L; RHO <- 0.25; PS <- c(500L, 1000L, 2000L)
ASTAR <- 0.30; BSTAR <- 0.15
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.10, 0.15, 0.20)
settings <- expand.grid(a = AGRID, b = BGRID)
setl <- data.frame(h = N^(-settings$b) / 2, alpha = N^(-settings$a))
IDX_PUB <- which(settings$a == ASTAR & settings$b == BSTAR)
MODEL_INDEX <- c(A1 = 1L, A2 = 3L)
OLD <- c(A1 = "M1", A2 = "M3")
KEEP <- 1:6   # campaign5 already used exactly these six rules

for (m in names(MODEL_INDEX)) for (p in PS) {
  cell_id <- match(p, c(2000L, 1000L, 500L))
  old <- readRDS(sprintf("results/campaign5/c5_n2000_p%04d_r025.rds", p))
  old <- Filter(function(z) z$model == OLD[[m]], old)
  for (r in c(1L, 7L, 42L)) {
    seed <- 131000021L + cell_id * 100003L + MODEL_INDEX[[m]] * 10007L +
      r * 307L
    d <- simulate_dataset5(N, p, RHO, m, seed)
    S <- matrix(NA_real_, p, nrow(setl)); R <- S
    for (k in seq_len(nrow(setl))) {
      S[, k] <- vapply(seq_len(p), function(j)
        score_coordinate_cpp(d$z[, j], d$y, setl$h[k], setl$alpha[k],
                             .05)[["score"]], numeric(1))
      R[, k] <- rank(S[, k], ties.method = "first")
    }
    amin <- apply(R, 1L, min)
    uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (N + 1))
    yu <- yu_score_matrix(uh, d$y, k = floor(0.072 * N), h = 1)
    ords <- c(list(order(S[, IDX_PUB], seq_len(p), na.last = TRUE),
                   order(amin, seq_len(p)),
                   order(-yu$scores, seq_len(p))),
              lapply(c(.90, .95, .99), function(tt)
                order(-qa_sis_scores(uh, d$y, tau = tt), seq_len(p))))
    rmax_new <- vapply(ords, function(o) max(match(A_GAMMA, o)), numeric(1))
    z <- old[[which(vapply(old, `[[`, 0L, "replicate") == r)]]
    stopifnot(all(rmax_new == z$rmax[KEEP]))
  }
  cat("PASS equivalence:", m, "== campaign5", OLD[[m]], "at p =", p, "\n")
}
cat("ALL EQUIVALENCE CHECKS PASS\n")
