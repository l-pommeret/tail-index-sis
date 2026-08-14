## Enforce a single canonical run for the design cell that is shared by the
## dimension, dependence and model tables: model A, n = p = 1000, rho = 0,
## a = b = 0.2. The dependence and models experiments already share one seed
## batch for that cell; the dimension experiment had drawn an independent
## batch, so the same cell was reported with two different Monte Carlo values.
## We keep the shared batch everywhere and rewrite the dimension summaries and
## the replicate-level metrics used by the figures.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: canonical_cell.R RESULTS_DIR FINAL_DIR")
rdir <- args[1L]; fdir <- args[2L]

cell <- function(z) z$model == "A" && z$n == 1000L && z$p == 1000L &&
  z$rho == 0 && z$a == 0.2 && z$b == 0.2

dep <- readRDS(file.path(rdir, "dependence.rds"))$jobs
canon <- Filter(cell, dep)
if (!length(canon)) stop("canonical cell not found in dependence.rds")
canon_seeds <- vapply(canon, `[[`, numeric(1), "seed")

summ <- function(w) {
  r <- vapply(w, function(q) max(q$active_ranks), numeric(1))
  data.frame(reps = length(w), sure20 = mean(r <= 20),
             exact = mean(vapply(w, function(q)
               max(q$active_ranks) <= length(q$active), logical(1))),
             mean_max_rank = mean(r),
             mean_under_rate = mean(unlist(lapply(w, `[[`, "under_rate"))),
             mean_elapsed = mean(vapply(w, `[[`, numeric(1), "elapsed")))
}

md_csv <- file.path(rdir, "tables", "main_dimension.csv")
md <- read.csv(md_csv)
i <- which(md$n == 1000 & md$p == 1000)
stopifnot(length(i) == 1L)
s <- summ(canon)
for (nm in names(s)) md[i, nm] <- s[[nm]]
write.csv(md[order(md$n, md$p), ], md_csv, row.names = FALSE)

d_grid <- c(4L, 10L, 20L, 50L, 100L, 200L, 500L)
byd_csv <- file.path(rdir, "tables", "main_dimension_by_d.csv")
byd <- read.csv(byd_csv)
for (dd in d_grid) {
  sel <- lapply(canon, function(z) z$ordering[seq_len(dd)])
  sure <- mean(mapply(function(z, s) all(z$active %in% s), canon, sel))
  mt <- mean(mapply(function(z, s) sum(z$active %in% s), canon, sel))
  ex <- mean(mapply(function(z, s) dd == length(z$active) &&
                      setequal(s, z$active), canon, sel))
  k <- which(byd$n == 1000 & byd$p == 1000 & byd$d == dd)
  if (length(k) == 1L) {
    byd$reps[k] <- length(canon); byd$sure[k] <- sure; byd$mean_true[k] <- mt
    byd$mean_max_rank[k] <- mean(vapply(canon, function(z)
      max(z$active_ranks), numeric(1)))
    byd$exact[k] <- ex
  }
}
write.csv(byd[order(byd$model, byd$n, byd$p, byd$rho, byd$d), ], byd_csv,
          row.names = FALSE)

## Replicate-level metrics feed the figures; drop the superseded batch so that
## the (n, p) = (1000, 1000) point is the canonical cell rather than a pool.
rm_csv <- file.path(fdir, "replicate_metrics.csv")
if (file.exists(rm_csv)) {
  raw <- read.csv(rm_csv)
  dup <- raw$model == "A" & raw$n == 1000 & raw$p == 1000 & raw$rho == 0 &
    raw$a == 0.2 & raw$b == 0.2 & !(raw$seed %in% canon_seeds)
  write.csv(raw[!dup, ], rm_csv, row.names = FALSE)
  cat("DROPPED", sum(dup), "superseded replicate rows\n")
}
cat("CANONICAL CELL sure20", s$sure20, "exact", s$exact,
    "mean_max_rank", s$mean_max_rank, "\n")
