## Tuning study for the Yoshida--Umezu competitor, mirroring the tuning
## study run for the proposed screen (code/R/tuning6.R).
##
## Why this exists.  In the previous comparison the competitor was run at
## the single tuning transplanted from its own paper (h = 1, k = 0.072 n),
## common to all four models, while the proposed screen was reported at a
## pair selected on this suite.  That is not a like-for-like comparison.
## Here the competitor is given a grid search of its own, and --- more
## generously still --- the best cell is kept SEPARATELY FOR EACH MODEL,
## whereas the proposed screen keeps one pair common to all four.
##
## The grid is chosen to be coherent with Yoshida--Umezu's own study:
##   * bandwidth h.  They work in a deliberate large-bandwidth regime,
##     report h in [0.04, 10], find h < 0.1 unusable, h = 1.5-2.0 best,
##     and h > 2 progressively oversmoothed.  The grid {0.5, 1, 1.5, 2, 3}
##     brackets their recommended region and both failure modes.
##   * intermediate sequence k.  Their baseline is k = 180 at n = 2500
##     (fraction 0.072), and their sensitivity study spans k in [121, 200],
##     i.e. fractions [0.048, 0.080].  The grid {0.05, 0.06, 0.072, 0.08}
##     covers that span at n = 2000.
## Their evaluation points and Epanechnikov kernel are kept unchanged.
##
## Common random numbers: the seed stream is EXACTLY tuning6's, so both
## methods are tuned on the same 200 datasets per model.
##
## usage: Rscript code/R/yu_tuning8.R OUTDIR [NREP] [CORES]
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "results/yu_tuning8"
NREP  <- if (length(args) >= 2) as.integer(args[2]) else 200L
CORES <- if (length(args) >= 3) as.integer(args[3]) else 7L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

source("code/R/generate5.R")
source("code/R/yoshida_umezu.R")
source("code/R/yoshida_umezu_grid.R")
source("code/R/yu_select.R")

N <- 2000L; P <- 1000L; RHO <- 0.25
HS <- c(0.5, 1, 1.5, 2, 3)
KFRAC <- c(0.05, 0.06, 0.072, 0.08)
KS <- as.integer(floor(KFRAC * N))
EVAL <- seq(.02, .98, length.out = 25)
MODELS <- c("A1", "A2", "A3", "B1")
BLOCK <- as.integer(Sys.getenv("YU8_BLOCK", "25"))

one_rep <- function(m, r, seed) {
  d <- simulate_dataset5(N, P, RHO, m, seed)
  u <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (N + 1))
  out <- vector("list", length(HS))
  for (ih in seq_along(HS)) {
    fit <- yu_score_matrix_grid(u, d$y, KS, HS[ih], EVAL)
    rmax <- vapply(seq_along(KS), function(ik) {
      o <- order(-fit$scores[, ik], seq_len(P), na.last = TRUE)
      max(match(A_GAMMA, o))
    }, numeric(1))
    out[[ih]] <- data.frame(model = m, replicate = r, seed = seed,
                            h = HS[ih], k = KS, k_fraction = KFRAC,
                            rmax = rmax,
                            undefined = colMeans(fit$undefined > 0))
  }
  do.call(rbind, out)
}

for (mi in seq_along(MODELS)) {
  m <- MODELS[mi]
  for (b0 in seq(1L, NREP, by = BLOCK)) {
    b1 <- min(NREP, b0 + BLOCK - 1L)
    path <- file.path(outdir, sprintf("yu8_%s_r%03d_%03d.rds", m, b0, b1))
    if (file.exists(path)) next
    out <- parallel::mclapply(b0:b1, function(r) {
      seed <- 141000041L + mi * 100019L + r * 211L   # tuning6 stream
      one_rep(m, r, seed)
    }, mc.cores = CORES, mc.preschedule = FALSE)
    bad <- vapply(out, function(z) inherits(z, "try-error") || is.null(z),
                  logical(1))
    if (any(bad)) stop("block ", path, ": ", sum(bad), " failed")
    saveRDS(out, path, compress = "xz")
    cat(format(Sys.time(), "%H:%M:%S"), m,
        sprintf("reps %d-%d done\n", b0, b1)); flush.console()
  }
}

## Summary and selection.
files <- list.files(outdir, pattern = "^yu8_", full.names = TRUE)
raw <- do.call(rbind, unlist(lapply(files, readRDS), recursive = FALSE))
tab <- do.call(rbind, lapply(
  split(raw, list(raw$model, raw$h, raw$k), drop = TRUE), function(z)
    data.frame(model = z$model[1], h = z$h[1], k = z$k[1],
               k_fraction = z$k_fraction[1], reps = nrow(z),
               sure4 = mean(z$rmax <= 4), sure20 = mean(z$rmax <= 20),
               ermax = round(mean(z$rmax), 3), medmax = median(z$rmax),
               undefined = mean(z$undefined))))
rownames(tab) <- NULL
tab <- tab[order(tab$model, tab$h, tab$k), ]
write.csv(tab, file.path(outdir, "yu_tuning_grid8.csv"), row.names = FALSE)

best <- yu_select_by_model(tab); rownames(best) <- NULL
write.csv(best, file.path(outdir, "yu_selected.csv"), row.names = FALSE)

cat("\nYU tuning grid (", nrow(tab), " cells,", NREP, "reps each)\n")
print(tab, row.names = FALSE)
cat("\nSelected per model (Sure-20 argmax, ties by E(Rmax)):\n")
print(best[, c("model", "h", "k", "k_fraction", "sure4", "sure20", "ermax")],
      row.names = FALSE)
cat("\nYU TUNING DONE\n")
