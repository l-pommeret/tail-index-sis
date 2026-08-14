args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: summarize_experiment.R INPUT.rds OUTPUT.csv")
x <- readRDS(args[1L])$jobs
key <- vapply(x, function(z) paste(z$model, z$signal, z$n, z$p, z$rho,
                                   z$a, z$b, sep = "|"), "")
groups <- split(x, key)
rows <- lapply(groups, function(w) {
  z <- w[[1L]]
  if (z$model == "N") {
    sure20 <- exact <- mean_max_rank <- NA_real_
  } else {
    sure20 <- mean(vapply(w, function(q) max(q$active_ranks) <= 20, logical(1)))
    exact <- mean(vapply(w, function(q) max(q$active_ranks) <= length(q$active), logical(1)))
    mean_max_rank <- mean(vapply(w, function(q) max(q$active_ranks), numeric(1)))
  }
  data.frame(model=z$model, signal=z$signal, n=z$n, p=z$p, rho=z$rho,
             a=z$a, b=z$b, reps=length(w), sure20=sure20, exact=exact,
             mean_max_rank=mean_max_rank,
             mean_under_rate=mean(unlist(lapply(w, `[[`, "under_rate"))),
             mean_elapsed=mean(vapply(w, `[[`, numeric(1), "elapsed")))
})
out <- do.call(rbind, rows)
dir.create(dirname(args[2L]), recursive = TRUE, showWarnings = FALSE)
write.csv(out[order(out$model, out$a, out$b),], args[2L], row.names = FALSE)
cat("WROTE", args[2L], "ROWS", nrow(out), "\n")
