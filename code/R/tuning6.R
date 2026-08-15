## Draft-5 tuning study, common-random-number design.
## For each model, generate NREP datasets once (seed stream below) and
## evaluate EVERY (a,b) cell of the grid on those same datasets, using the
## multi-setting streaming generator (one data pass, 42 score passes).
## Grid: a in {.25,...,.50} x b in {0,...,.40}; alpha = n^-a, h = n^-b / 2.
## Seed stream: 141000041 + model_index*100019 + r*211, disjoint from all
## previous streams.  Checkpointed per (model, replication-block).
## usage: Rscript code/R/tuning6.R OUTDIR [NREP] [CORES]
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "results/tuning6"
NREP  <- if (length(args) >= 2) as.integer(args[2]) else 200L
CORES <- if (length(args) >= 3) as.integer(args[3]) else 6L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate5.R")

N <- 2000L; P <- 1000L; RHO <- 0.25
AS <- c(.25, .30, .35, .40, .45, .50)
BS <- c(0, .05, .10, .15, .20, .30, .40)
grid <- expand.grid(a = AS, b = BS)
settings <- data.frame(h = N^(-grid$b) / 2, alpha = N^(-grid$a))

MODELS <- c("A1", "A2", "A3", "B1")
BLOCK <- 25L  # replications per checkpoint file
for (mi in seq_along(MODELS)) {
  m <- MODELS[mi]
  for (b0 in seq(1L, NREP, by = BLOCK)) {
    b1 <- min(NREP, b0 + BLOCK - 1L)
    path <- file.path(outdir, sprintf("t6_%s_r%03d_%03d.rds", m, b0, b1))
    if (file.exists(path)) next
    out <- parallel::mclapply(b0:b1, function(r) {
      seed <- 141000041L + mi * 100019L + r * 211L
      S <- simulate_scores_streaming5(N, P, RHO, m, seed, settings)
      rmax <- vapply(seq_len(nrow(grid)), function(k) {
        o <- order(S[, k], seq_len(P), na.last = TRUE)
        max(match(A_GAMMA, o))
      }, numeric(1))
      list(model = m, replicate = r, seed = seed, rmax = rmax)
    }, mc.cores = CORES, mc.preschedule = FALSE)
    bad <- vapply(out, function(z) inherits(z, "try-error") || is.null(z),
                  logical(1))
    if (any(bad)) stop("block ", path, ": ", sum(bad), " failed")
    saveRDS(out, path, compress = "xz")
    cat(format(Sys.time(), "%H:%M:%S"), m, sprintf("reps %d-%d done\n", b0, b1))
    flush.console()
  }
}

## Summary: Sure-20 per (model, a, b).
files <- list.files(outdir, pattern = "^t6_", full.names = TRUE)
jobs <- unlist(lapply(files, readRDS), recursive = FALSE)
res <- do.call(rbind, lapply(jobs, function(z)
  data.frame(model = z$model, a = grid$a, b = grid$b, rmax = z$rmax)))
ag <- aggregate(rmax <= 20 ~ model + a + b, res, mean)
names(ag)[4] <- "sure20"
write.csv(ag, file.path(outdir, "summary.csv"), row.names = FALSE)
w <- reshape(ag, idvar = c("a", "b"), timevar = "model", direction = "wide")
w$mean4 <- rowMeans(w[, grep("sure20", names(w))])
print(w[order(-w$mean4), ][1:10, ], digits = 3)
cat("TUNING6 DONE\n")
