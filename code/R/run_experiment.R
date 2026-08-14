args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: run_experiment.R CONFIG.csv OUTPUT.rds [CORES]")
config_path <- args[1L]
output_path <- args[2L]
cores <- if (length(args) >= 3L) as.integer(args[3L]) else 1L

source("code/R/generate.R")
source("code/R/streaming.R")
library(Rcpp)
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)

cfg <- read.csv(config_path, stringsAsFactors = FALSE)
jobs <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(i) {
  z <- cfg[rep(i, cfg$reps[i]), , drop = FALSE]
  z$replicate <- seq_len(cfg$reps[i])
  z$config_row <- i
  z
}))
jobs$seed <- 1000003L + jobs$config_row * 10007L + jobs$replicate * 101L

one_job <- function(z) {
  started <- proc.time()[[3L]]
  h <- z$n^(-z$b) / 2 # q=1 and uniform-norm unit-ball volume 2
  alpha <- z$n^(-z$a)
  d <- simulate_score_streaming(z$n, z$p, z$rho, z$model, z$signal, z$seed,
                                h, alpha, z$epsilon, z$block_size)
  scores <- d$scores
  ordering <- order(scores, seq_along(scores), na.last = TRUE)
  active <- if (z$model == "N") integer() else if (z$model == "E") 1:2 else 1:4
  active_ranks <- if (length(active)) match(active, ordering) else integer()
  elapsed <- proc.time()[[3L]] - started
  list(config_row = z$config_row, replicate = z$replicate, seed = z$seed,
       model = z$model, signal = z$signal, n = z$n, p = z$p, rho = z$rho,
       a = z$a, b = z$b, alpha = alpha, h = h, epsilon = z$epsilon,
       scores = scores, ordering = ordering, active = active,
       active_ranks = active_ranks,
       under_rate = d$under_rate, mean_local_n = d$mean_local_n,
       mean_local_k = d$mean_local_k, y_min = d$y_min,
       y_finite = d$y_finite, gamma_range = d$gamma_range, elapsed = elapsed)
}

if (.Platform$OS.type == "unix" && cores > 1L) {
  out <- parallel::mclapply(seq_len(nrow(jobs)), function(i) one_job(jobs[i,]),
                            mc.cores = cores, mc.preschedule = FALSE)
} else out <- lapply(seq_len(nrow(jobs)), function(i) one_job(jobs[i,]))

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
meta <- list(config_path = normalizePath(config_path), config = cfg,
             started_utc = format(Sys.time(), tz = "UTC"),
             r_version = R.version.string, rcpp_version = as.character(packageVersion("Rcpp")),
             cores = cores)
saveRDS(list(meta = meta, jobs = out), output_path, compress = "xz")
cat("WROTE", output_path, "JOBS", length(out), "\n")
