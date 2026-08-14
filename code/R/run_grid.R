## Extended method-comparison grid for the Draft-3 designs.
##
## Same estimators, tuning and data-generating processes as
## `code/R/run_draft3.R comparison`, but over an enlarged design space:
##   models  M1, M2, M3, M4
##   n       1000, 2000, 5000
##   p       200, 500, 1000, 2000
##   rho     0, 0.20, 0.25, 0.30, 0.40, 0.50
## with REPS Monte Carlo replications per cell (default 40).
##
## Seed stream: 12000019 + cell*10007 + r*101, disjoint from the pilot
## (7000003+...), tuning-M2 (2200003+...), comparison (4000037+...),
## comparator (7000003+i*10007+...) and real-data (5e5, 9e5) streams.
##
## Checkpointed per cell and resumable; cells are scheduled longest-first and
## dispatched replicate-by-replicate over `cores` workers.
##
## usage: KAPPA=0.20 Rscript code/R/run_grid.R OUTDIR CORES [REPS]

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[1L] else "results/grid/comparison_cells"
cores  <- if (length(args) >= 2L) as.integer(args[2L]) else 90L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R")
source("code/R/qa_sis.R")
source("code/R/yoshida_umezu.R")

EPS <- 0.05
stopifnot(nchar(Sys.getenv("KAPPA")) > 0)
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
TAUS  <- c(0.90, 0.95, 0.975, 0.99)
REPS  <- if (length(args) >= 3L) as.integer(args[3L]) else 40L

## The four axes can be narrowed with env vars (comma separated), for smoke
## tests and partial reruns; the defaults are the full design space.
axis <- function(var, default) {
  v <- Sys.getenv(var)
  if (!nzchar(v)) return(default)
  out <- strsplit(v, ",", fixed = TRUE)[[1]]
  if (is.numeric(default)) storage.mode(out) <- storage.mode(default)
  out
}
MODELS <- axis("GRID_MODELS", c("M1", "M2", "M3", "M4"))
NS     <- axis("GRID_N",   c(1000L, 2000L, 5000L))
PS     <- axis("GRID_P",   c(200L, 500L, 1000L, 2000L))
RHOS   <- axis("GRID_RHO", c(0, 0.20, 0.25, 0.30, 0.40, 0.50))

grid <- expand.grid(rho = RHOS, p = PS, n = NS, model = MODELS,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid$cell <- seq_len(nrow(grid))
grid$id <- sprintf("g_%s_n%04d_p%04d_r%03.0f", grid$model, grid$n, grid$p,
                   grid$rho * 100)
## per-replicate cost proxy in seconds (code/R/bench_grid.R), for scheduling
cost1 <- c("1000" = 8.0, "2000" = 16.7, "5000" = 68.7)
grid$cost <- REPS * cost1[as.character(grid$n)] * grid$p / 1000

one_rep <- function(g, r) {
  seed <- 12000019L + g$cell * 10007L + r * 101L
  n <- g$n; p <- g$p
  d <- simulate_dataset3(n, p, g$rho, g$model, seed)
  d$u <- NULL; d$gamma <- NULL          # population PIT never reaches a screen
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (n + 1))
  h <- n^(-BSTAR) / 2; alpha <- n^(-ASTAR)
  t0 <- proc.time()[3]
  sc <- vapply(seq_len(p), function(j)
    score_coordinate_cpp(d$z[, j], d$y, h, alpha, EPS), numeric(5))["score", ]
  t_ours <- proc.time()[3] - t0
  o_ours <- order(sc, seq_len(p), na.last = TRUE)
  t0 <- proc.time()[3]
  yu <- yu_score_matrix(uh, d$y, k = floor(0.072 * n), h = 1)
  t_yu <- proc.time()[3] - t0
  o_yu <- order(-yu$scores, seq_len(p))
  t0 <- proc.time()[3]
  qs <- lapply(TAUS, function(tt) {
    o <- order(-qa_sis_scores(uh, d$y, tau = tt), seq_len(p))
    list(ranks = match(1:4, o), top4 = o[1:4])
  })
  t_q <- proc.time()[3] - t0
  names(qs) <- paste0("q", TAUS * 1000)
  list(model = g$model, n = n, p = p, rho = g$rho, replicate = r, seed = seed,
       ranks_ours = match(1:4, o_ours), top4_ours = o_ours[1:4],
       ranks_yu = match(1:4, o_yu), top4_yu = o_yu[1:4],
       yu_undefined = mean(!is.finite(yu$scores)), q = qs,
       elapsed = c(ours = unname(t_ours), yu = unname(t_yu),
                   qsis = unname(t_q) / length(TAUS)))
}

run_batch <- function(rows) {
  jobs <- do.call(rbind, lapply(seq_len(nrow(rows)), function(k)
    data.frame(k = k, r = seq_len(REPS))))
  out <- parallel::mclapply(seq_len(nrow(jobs)), function(ix)
    try(one_rep(rows[jobs$k[ix], ], jobs$r[ix]), silent = TRUE),
    mc.cores = cores, mc.preschedule = FALSE)
  for (k in seq_len(nrow(rows))) {
    sel <- out[jobs$k == k]
    bad <- vapply(sel, function(z) inherits(z, "try-error") || is.null(z),
                  logical(1))
    if (any(bad)) {
      msg <- unlist(lapply(sel[bad], as.character))[1]
      stop("cell ", rows$id[k], ": ", sum(bad), " replicates failed: ", msg)
    }
    saveRDS(sel, file.path(outdir, paste0(rows$id[k], ".rds")), compress = "xz")
    cat(format(Sys.time(), "%H:%M:%S"), "cell", rows$id[k], "done\n")
  }
  flush.console()
}

done <- file.exists(file.path(outdir, paste0(grid$id, ".rds")))
todo <- grid[!done, ]
todo <- todo[order(-todo$cost), ]
cat(sprintf("%d cells total, %d already on disk, %d to run (%.1f core-hours)\n",
            nrow(grid), sum(done), nrow(todo), sum(todo$cost) / 3600))
flush.console()

batch_budget <- cores * 420          # aim for ~7 minutes of wall time per batch
start <- proc.time()[3]
while (nrow(todo) > 0L) {
  take <- max(1L, sum(cumsum(todo$cost) <= batch_budget))
  rows <- todo[seq_len(take), ]
  cat(sprintf("%s batch of %d cells (%.2f core-hours), %d cells left after\n",
              format(Sys.time(), "%H:%M:%S"), nrow(rows), sum(rows$cost) / 3600,
              nrow(todo) - take))
  flush.console()
  run_batch(rows)
  todo <- todo[-seq_len(take), ]
  cat(sprintf("%s elapsed %.1f min, remaining budget %.1f core-hours\n",
              format(Sys.time(), "%H:%M:%S"), (proc.time()[3] - start) / 60,
              sum(todo$cost) / 3600))
  flush.console()
}
cat("ALL CELLS DONE\n")
