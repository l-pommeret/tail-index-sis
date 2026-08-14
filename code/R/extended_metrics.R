args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: extended_metrics.R INPUT.rds OUTPUT.csv")
x <- readRDS(args[1L])$jobs
d_grid <- c(4L, 10L, 20L, 50L, 100L, 200L, 500L)
rows <- list(); at <- 1L
for (z in x) {
  if (!length(z$active)) next
  for (dd in unique(pmin(d_grid, z$p))) {
    selected <- z$ordering[seq_len(dd)]
    rows[[at]] <- data.frame(model=z$model, signal=z$signal, n=z$n, p=z$p,
    rho=z$rho, a=z$a, b=z$b, replicate=z$replicate, d=dd,
      sure=all(z$active %in% selected), true_in_top=sum(z$active %in% selected),
      max_active_rank=max(z$active_ranks), exact=(dd==length(z$active) &&
        setequal(selected,z$active)))
    at <- at + 1L
  }
}
raw <- do.call(rbind, rows)
key <- interaction(raw[,c("model","signal","n","p","rho","a","b","d")], drop=TRUE)
out <- do.call(rbind, lapply(split(raw,key), function(w) data.frame(
  model=w$model[1], signal=w$signal[1], n=w$n[1], p=w$p[1], rho=w$rho[1],
  a=w$a[1], b=w$b[1], d=w$d[1], reps=nrow(w),
  sure=mean(w$sure), mean_true=mean(w$true_in_top),
  mean_max_rank=mean(w$max_active_rank), exact=mean(w$exact))))
dir.create(dirname(args[2L]), recursive=TRUE, showWarnings=FALSE)
write.csv(out[order(out$model,out$n,out$p,out$rho,out$d),],args[2L],row.names=FALSE)
