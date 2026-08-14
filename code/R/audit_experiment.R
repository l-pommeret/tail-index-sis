args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: audit_experiment.R OUTPUT.rds")
x <- readRDS(args[1L]); jobs <- x$jobs
seeds <- vapply(jobs, `[[`, integer(1), "seed")
if (isTRUE(x$meta$common_random_numbers)) {
  cfg_rows <- vapply(jobs, `[[`, integer(1), "config_row")
  seed_sets <- split(seeds, cfg_rows)
  stopifnot(all(vapply(seed_sets, identical, logical(1), seed_sets[[1L]])),
            length(unique(seeds)) == unique(x$meta$config$reps))
} else stopifnot(length(unique(seeds)) == length(jobs))
stopifnot(length(jobs) == sum(x$meta$config$reps),
          all(vapply(jobs, `[[`, logical(1), "y_finite")),
          all(vapply(jobs, `[[`, numeric(1), "y_min") > 1),
          all(vapply(jobs, function(z) length(z$scores) == z$p, logical(1))),
          all(vapply(jobs, function(z) all(is.finite(z$scores)), logical(1))),
          all(vapply(jobs, function(z) identical(sort(z$ordering), seq_len(z$p)), logical(1))))
cat("PASS", args[1L], "jobs", length(jobs), "\n")
