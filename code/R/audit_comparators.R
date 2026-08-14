args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: audit_comparators.R OUTPUT.rds")
x <- readRDS(args[1L]); jobs <- x$jobs
stopifnot(length(jobs) == sum(x$meta$config$reps),
          length(unique(vapply(jobs, `[[`, integer(1), "seed"))) == length(jobs),
          all(vapply(jobs, function(z)
            length(z$ours_scores) == z$p && length(z$yu_scores) == z$p &&
            all(is.finite(z$ours_scores)) && all(is.finite(z$yu_scores)) &&
            identical(sort(z$ours_ordering), seq_len(z$p)) &&
            identical(sort(z$yu_ordering), seq_len(z$p)) &&
            isTRUE(z$yu_empirical_pit), logical(1))))
cat("PASS", args[1L], "jobs", length(jobs), "finite scores",
    2L * sum(vapply(jobs, `[[`, integer(1), "p")),
    "YU_EMPIRICAL_PIT TRUE\n")
