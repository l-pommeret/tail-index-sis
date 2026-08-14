## Draft-3 experiment driver (checkpointed, resumable).
## KAPPA env var must be set to the pilot-selected value.
## Seed streams (all disjoint from the pilot stream 7000003+...):
##   tuning M2 (draft-3): 2200003 + i*100003 + r*211
##   comparison:          4000037 + i*100043 + r*211  (M1/M3/M4 datasets
##                        identical to Draft 2's comparison cells)
##   YU bandwidth grid:   same comparison seeds (paired).
## usage: Rscript code/R/run_draft3.R {tuning_m2|comparison|yugrid_m2} OUTDIR CORES
args <- commandArgs(trailingOnly = TRUE)
what <- args[1L]; outdir <- args[2L]; cores <- as.integer(args[3L])
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R")
EPS <- 0.05
stopifnot(nchar(Sys.getenv("KAPPA")) > 0)

run_cell <- function(cell_id, reps, fn) {
  path <- file.path(outdir, paste0(cell_id, ".rds"))
  if (file.exists(path)) return(invisible(NULL))
  out <- parallel::mclapply(seq_len(reps), fn, mc.cores = cores,
                            mc.preschedule = FALSE)
  bad <- vapply(out, function(z) inherits(z, "try-error") || is.null(z),
                logical(1))
  if (any(bad)) stop("cell ", cell_id, ": ", sum(bad), " failed")
  saveRDS(out, path, compress = "xz")
  cat(format(Sys.time(), "%H:%M:%S"), "cell", cell_id, "done\n"); flush.console()
}

if (what == "tuning_m2") {
  grid <- expand.grid(a = c(.25,.30,.35,.40,.45,.50),
                      b = c(0,.05,.10,.15,.20,.30,.40))
  grid <- grid[order(grid$b), ]
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    cid <- sprintf("tun3_M2_a%02.0f_b%02.0f", g$a*100, g$b*100)
    n <- 2000L; h <- n^(-g$b)/2; alpha <- n^(-g$a)
    seed0 <- 2200003L + i * 100003L
    run_cell(cid, 200L, function(r) {
      seed <- seed0 + r * 211L
      sc <- simulate_score_streaming3(n, 1000L, 0.25, "M2", seed, h, alpha, EPS)
      o <- order(sc, seq_along(sc), na.last = TRUE)
      list(a = g$a, b = g$b, replicate = r, seed = seed,
           active_ranks = match(1:4, o))
    })
  }
} else if (what == "comparison") {
  source("code/R/qa_sis.R"); source("code/R/yoshida_umezu.R")
  ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
  BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
  TAUS <- c(0.90, 0.95, 0.975, 0.99)
  for (i in 1:4) {
    m <- c("M1","M2","M3","M4")[i]
    cid <- paste0("cmp3_", m)
    seed0 <- 4000037L + i * 100043L
    n <- 2000L; p <- 1000L
    h <- n^(-BSTAR)/2; alpha <- n^(-ASTAR)
    run_cell(cid, 200L, function(r) {
      seed <- seed0 + r * 211L
      d <- simulate_dataset3(n, p, 0.25, m, seed)
      uh <- apply(d$z, 2L, function(x) rank(x, ties.method="average")/(n+1))
      sc <- vapply(seq_len(p), function(j)
        score_coordinate_cpp(d$z[,j], d$y, h, alpha, EPS), numeric(5))["score",]
      o_ours <- order(sc, seq_len(p), na.last = TRUE)
      yu <- yu_score_matrix(uh, d$y, k = floor(0.072*n), h = 1)
      o_yu <- order(-yu$scores, seq_len(p))
      qs <- lapply(TAUS, function(tt) {
        o <- order(-qa_sis_scores(uh, d$y, tau = tt), seq_len(p))
        list(ranks = match(1:4, o), top4 = o[1:4])
      })
      names(qs) <- paste0("q", TAUS*1000)
      list(model = m, replicate = r, seed = seed,
           ranks_ours = match(1:4, o_ours), top4_ours = o_ours[1:4],
           ranks_yu = match(1:4, o_yu), top4_yu = o_yu[1:4],
           q = qs)
    })
  }
} else if (what == "yugrid_m2") {
  source("code/R/yoshida_umezu.R")
  seed0 <- 4000037L + 2L * 100043L
  n <- 2000L; p <- 1000L
  run_cell("yu3_M2", 200L, function(r) {
    seed <- seed0 + r * 211L
    d <- simulate_dataset3(n, p, 0.25, "M2", seed)
    uh <- apply(d$z, 2L, function(x) rank(x, ties.method="average")/(n+1))
    res <- lapply(c(0.5, 1.5), function(hh) {
      yu <- yu_score_matrix(uh, d$y, k = floor(0.072*n), h = hh)
      match(1:4, order(-yu$scores, seq_len(p)))
    })
    list(replicate = r, seed = seed,
         ranks_h05 = res[[1]], ranks_h15 = res[[2]])
  })
} else stop("unknown")
cat("ALL CELLS DONE:", what, "\n")
