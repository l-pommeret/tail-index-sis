## Draft-2 driver with per-cell checkpointing (kill-safe, resumable).
## usage: Rscript code/R/run_draft2b.R {tuning|sensitivity|comparison} OUTDIR CORES
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
what <- args[1L]; outdir <- args[2L]; cores <- as.integer(args[3L])
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate2.R")
EPS <- 0.05; BLOCK <- 256L
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.15"))

run_cell <- function(cell_id, reps, seed0, fn) {
  path <- file.path(outdir, paste0(cell_id, ".rds"))
  if (file.exists(path)) return(invisible(NULL))
  out <- parallel::mclapply(seq_len(reps), function(r) fn(r, seed0 + r * 211L),
                            mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(out, function(z) inherits(z, "try-error") || is.null(z), logical(1))
  if (any(bad)) stop("cell ", cell_id, ": ", sum(bad), " failed")
  saveRDS(out, path, compress = "xz")
  cat(format(Sys.time(), "%H:%M:%S"), "cell", cell_id, "done\n"); flush.console()
}

if (what == "tuning") {
  grid <- expand.grid(model = c("M1","M2","M3","M4"),
                      a = c(.25,.30,.35,.40,.45,.50),
                      b = c(0,.05,.10,.15,.20,.30,.40),
                      stringsAsFactors = FALSE)
  ## cheap cells first so most of the grid lands early
  grid <- grid[order(grid$b), ]
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    cid <- sprintf("tun_%s_a%02.0f_b%02.0f", g$model, g$a*100, g$b*100)
    n <- 2000L; h <- n^(-g$b)/2; alpha <- n^(-g$a)
    seed0 <- 2000003L + i * 100003L
    run_cell(cid, 200L, seed0, function(r, seed) {
      d <- simulate_score_streaming2(n, 1000L, 0.25, g$model, seed, h, alpha, EPS, BLOCK)
      o <- order(d$scores, seq_along(d$scores), na.last = TRUE)
      list(model = g$model, a = g$a, b = g$b, replicate = r, seed = seed,
           active_ranks = match(1:4, o))
    })
  }
} else if (what == "sensitivity") {
  cells <- rbind(
    data.frame(n = c(1000L, 2000L, 5000L), p = 1000L, rho = 0.25),
    data.frame(n = 2000L, p = c(100L, 500L, 5000L), rho = 0.25),
    data.frame(n = 2000L, p = 1000L, rho = c(0, 0.5, 0.75)))
  for (i in seq_len(nrow(cells))) {
    g <- cells[i, ]
    cid <- sprintf("sen_n%d_p%d_r%03.0f", g$n, g$p, g$rho*100)
    h <- g$n^(-BSTAR)/2; alpha <- g$n^(-ASTAR)
    seed0 <- 3000017L + i * 100019L
    run_cell(cid, 200L, seed0, function(r, seed) {
      d <- simulate_score_streaming2(g$n, g$p, g$rho, "M1", seed, h, alpha, EPS, BLOCK)
      o <- order(d$scores, seq_along(d$scores), na.last = TRUE)
      list(n = g$n, p = g$p, rho = g$rho, replicate = r, seed = seed,
           a = ASTAR, b = BSTAR, active_ranks = match(1:4, o))
    })
  }
} else if (what == "comparison") {
  source("code/R/qa_sis.R"); source("code/R/yoshida_umezu.R")
  for (i in 1:4) {
    m <- c("M1","M2","M3","M4")[i]
    cid <- paste0("cmp_", m)
    seed0 <- 4000037L + i * 100043L
    n <- 2000L; p <- 1000L
    h <- n^(-BSTAR)/2; alpha <- n^(-ASTAR)
    run_cell(cid, 200L, seed0, function(r, seed) {
      d <- simulate_dataset2(n, p, 0.25, m, seed, BLOCK)
      uh <- apply(d$z, 2L, function(x) rank(x, ties.method="average")/(n+1))
      sc <- vapply(seq_len(p), function(j)
        score_coordinate_cpp(d$z[,j], d$y, h, alpha, EPS), numeric(5))["score",]
      o_ours <- order(sc, seq_len(p), na.last = TRUE)
      yu <- yu_score_matrix(uh, d$y, k = floor(0.072*n), h = 1)
      o_yu <- order(-yu$scores, seq_len(p))
      q95 <- qa_sis_scores(uh, d$y, tau = 0.95); o_q95 <- order(-q95, seq_len(p))
      q99 <- qa_sis_scores(uh, d$y, tau = 0.99); o_q99 <- order(-q99, seq_len(p))
      rk <- function(o) match(1:4, o)
      list(model = m, replicate = r, seed = seed,
           ranks_ours = rk(o_ours), ranks_yu = rk(o_yu),
           ranks_q95 = rk(o_q95), ranks_q99 = rk(o_q99),
           top20_q95 = o_q95[1:20], top20_q99 = o_q99[1:20])
    })
  }
} else stop("unknown experiment")
cat("ALL CELLS DONE:", what, "\n")
