## Emit the pilot table (all tested kappa) for the supplement.
res <- readRDS("results/draft3/pilot_m2.rds")
rows <- lapply(names(res), function(kap) {
  x <- res[[kap]]
  rt <- sapply(x, function(z) max(z$rt)); rq <- sapply(x, function(z) max(z$rq))
  data.frame(kappa = as.numeric(kap), reps = length(x),
    tail_sure4 = mean(rt <= 4), tail_sure20 = mean(rt <= 20),
    tail_ermax = round(mean(rt), 1),
    q95_sure4 = mean(rq <= 4), q95_sure20 = mean(rq <= 20),
    q95_ermax = round(mean(rq), 1),
    tail_top4_gamma = mean(sapply(x, function(z) sum(z$top4_t %in% 1:4))),
    tail_top4_scale = mean(sapply(x, function(z) sum(z$top4_t %in% 5:8))),
    q95_top4_gamma = mean(sapply(x, function(z) sum(z$top4_q %in% 1:4))),
    q95_top4_scale = mean(sapply(x, function(z) sum(z$top4_q %in% 5:8))))
})
out <- do.call(rbind, rows); out <- out[order(out$kappa), ]
write.csv(out, "results/draft3/pilot_summary.csv", row.names = FALSE)
print(out, digits = 3)
