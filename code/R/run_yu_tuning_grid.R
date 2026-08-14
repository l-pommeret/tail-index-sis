args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else "results/yu_tuning_grid.rds"
cores <- if (length(args) >= 2L) as.integer(args[2L]) else 1L
source("code/R/generate.R")
source("code/R/yoshida_umezu.R")
cfg <- read.csv("config/comparators.csv", stringsAsFactors = FALSE)
tuning <- expand.grid(k_fraction = c(.05, .072, .10), h = c(.5, 1, 1.5))
jobs <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(i)
  data.frame(config_row = i, replicate = seq_len(cfg$reps[i]))))
one <- function(ii) {
  row <- jobs[ii, ]; z <- cfg[row$config_row, ]
  seed <- 7000003L + row$config_row * 10007L + row$replicate * 101L
  d <- simulate_dataset(z$n, z$p, z$rho, z$model, z$signal, seed, z$block_size)
  u <- rank_uniform_matrix(d$x)
  active <- if (z$model == "E") 1:2 else 1:4
  do.call(rbind, lapply(seq_len(nrow(tuning)), function(g) {
    k <- floor(tuning$k_fraction[g] * z$n)
    fit <- yu_score_matrix(u, d$y, k, tuning$h[g], seq(.02, .98, length.out = 25))
    ord <- order(-fit$scores, seq_along(fit$scores), na.last = TRUE)
    ranks <- match(active, ord)
    data.frame(config_row = row$config_row, replicate = row$replicate, seed = seed,
      model = z$model, rho = z$rho, k_fraction = tuning$k_fraction[g], h = tuning$h[g],
      sure20 = max(ranks) <= 20, exact = max(ranks) <= length(active),
      mean_max_rank = max(ranks), undefined = mean(!is.finite(fit$scores)))
  }))
}
ans <- if (.Platform$OS.type == "unix" && cores > 1L)
  parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = cores, mc.preschedule = FALSE) else
  lapply(seq_len(nrow(jobs)), one)
raw <- do.call(rbind, ans)
summary <- aggregate(cbind(sure20, exact, mean_max_rank, undefined) ~
  model + rho + k_fraction + h, raw, mean)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(config = cfg, tuning = tuning, raw = raw, summary = summary), out, compress = "xz")
write.csv(summary, sub("\\.rds$", ".csv", out), row.names = FALSE)
cat("PASS YU tuning grid", nrow(raw), "fits\n")
