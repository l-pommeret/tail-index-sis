args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L)
  stop("usage: run_alpha_common_random.R CONFIG.csv OUTPUT.rds [CORES]")
cores <- if (length(args) >= 3L) as.integer(args[3L]) else 1L

source("code/R/generate.R")
source("code/R/streaming.R")
library(Rcpp)
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)

cfg <- read.csv(args[1L], stringsAsFactors = FALSE)
stopifnot(nrow(cfg) > 1L, length(unique(cfg$reps)) == 1L,
          length(unique(cfg$model)) == 1L, length(unique(cfg$signal)) == 1L,
          length(unique(cfg$n)) == 1L, length(unique(cfg$p)) == 1L,
          length(unique(cfg$rho)) == 1L, length(unique(cfg$b)) == 1L,
          length(unique(cfg$epsilon)) == 1L,
          length(unique(cfg$block_size)) == 1L)

jobs <- expand.grid(config_row = seq_len(nrow(cfg)),
                    replicate = seq_len(cfg$reps[1L]))
# Common random numbers: seed depends only on replicate, never on alpha row.
jobs$seed <- 9000001L + jobs$replicate * 101L

one_job <- function(i) {
  jj <- jobs[i, ]; z <- cfg[jj$config_row, ]
  started <- proc.time()[[3L]]
  alpha <- z$n^(-z$a); h <- z$n^(-z$b) / 2
  d <- simulate_score_streaming(z$n, z$p, z$rho, z$model, z$signal,
                                jj$seed, h, alpha, z$epsilon, z$block_size)
  ordering <- order(d$scores, seq_along(d$scores), na.last = TRUE)
  active <- 1:4
  list(config_row = jj$config_row, replicate = jj$replicate, seed = jj$seed,
       model = z$model, signal = z$signal, n = z$n, p = z$p, rho = z$rho,
       a = z$a, b = z$b, alpha = alpha, h = h, epsilon = z$epsilon,
       scores = d$scores, ordering = ordering, active = active,
       active_ranks = match(active, ordering), under_rate = d$under_rate,
       mean_local_n = d$mean_local_n, mean_local_k = d$mean_local_k,
       y_min = d$y_min, y_finite = d$y_finite,
       gamma_range = d$gamma_range,
       elapsed = proc.time()[[3L]] - started)
}

if (.Platform$OS.type == "unix" && cores > 1L) {
  out <- parallel::mclapply(seq_len(nrow(jobs)), one_job, mc.cores = cores,
                            mc.preschedule = FALSE)
} else {
  out <- lapply(seq_len(nrow(jobs)), one_job)
}

# Verify the common-random-number design at archive creation time.
seed_sets <- split(vapply(out, `[[`, integer(1), "seed"),
                   vapply(out, `[[`, integer(1), "config_row"))
stopifnot(all(vapply(seed_sets, identical, logical(1), seed_sets[[1L]])))
meta <- list(config_path = normalizePath(args[1L]), config = cfg,
             common_random_numbers = TRUE,
             seed_formula = "9000001 + replicate * 101",
             started_utc = format(Sys.time(), tz = "UTC"),
             r_version = R.version.string,
             rcpp_version = as.character(packageVersion("Rcpp")), cores = cores)
saveRDS(list(meta = meta, jobs = out), args[2L], compress = "xz")
cat("WROTE", args[2L], "JOBS", length(out), "COMMON_RANDOM_NUMBERS TRUE\n")
