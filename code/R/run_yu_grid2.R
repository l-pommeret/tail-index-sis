## Fairness supplement: Yoshida-Umezu at h in {0.5, 1, 1.5} (k = 0.072 n fixed)
## on the SAME comparison replications: seed increment 211 reproduces the
## exact datasets of cmp_*.rds.
## usage: Rscript code/R/run_yu_grid2.R OUTDIR CORES
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1L]; cores <- as.integer(args[2L])
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
source("code/R/generate2.R"); source("code/R/yoshida_umezu.R")
for (i in 1:4) {
  m <- c("M1","M2","M3","M4")[i]
  path <- file.path(outdir, paste0("yu_", m, ".rds"))
  if (file.exists(path)) next
  seed0 <- 4000037L + i * 100043L
  out <- parallel::mclapply(1:200, function(r) {
    seed <- seed0 + r * 211L  # matches cmp_*.rds datasets exactly
    n <- 2000L; p <- 1000L
    d <- simulate_dataset2(n, p, 0.25, m, seed, 256L)
    uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average")/(n+1))
    res <- lapply(c(0.5, 1.5), function(hh) {
      yu <- yu_score_matrix(uh, d$y, k = floor(0.072*n), h = hh)
      match(1:4, order(-yu$scores, seq_len(p)))
    })
    list(model = m, replicate = r, seed = seed,
         ranks_h05 = res[[1]], ranks_h15 = res[[2]])
  }, mc.cores = cores, mc.preschedule = FALSE)
  saveRDS(out, path, compress = "xz")
  cat(format(Sys.time(), "%H:%M:%S"), "yu grid", m, "done\n")
}
cat("YU GRID DONE\n")
