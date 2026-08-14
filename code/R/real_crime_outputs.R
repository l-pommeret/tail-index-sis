## Draft-2 real-data figures and summaries from results/draft2/real_crime.rds.
x <- readRDS("results/draft2/real_crime.rds")
p <- x$p; nm <- x$var_names

primary_rank <- rank(x$primary_scores, ties.method = "first")
lead <- order(x$primary_scores)[1:10]

## Rank stability across the 3x3 tuning grid.
grid_mat <- vapply(x$grid_ranks, function(r) r[lead], numeric(10))
## Half-sample rank distribution for the leaders.
half_mat <- vapply(x$half_ranks, function(r) r[lead], numeric(10))
## Tie-audit ranks.
tie_mat <- vapply(x$tie_ranks, function(r) r[lead], numeric(10))

summ <- data.frame(
  variable = nm[lead],
  primary = primary_rank[lead],
  grid_med = apply(grid_mat, 1, quantile, .5, type = 1),
  grid_max = apply(grid_mat, 1, max),
  half_med = apply(half_mat, 1, quantile, .5, type = 1),
  half_q10 = apply(half_mat, 1, quantile, .1, type = 1),
  half_q90 = apply(half_mat, 1, quantile, .9, type = 1),
  top5_freq = rowMeans(half_mat <= 5),
  tie_max = apply(tie_mat, 1, max),
  yu_rank = rank(-x$yu_scores, ties.method = "first")[lead],
  q95_rank = rank(-x$q95_scores, ties.method = "first")[lead])
write.csv(summ, "results/draft2/real_crime_leaders.csv", row.names = FALSE)
print(summ, digits = 3)
cat("\nn =", x$n, " p =", p, " alpha =", round(x$alpha, 4),
    " h =", round(x$h, 4), " n*alpha =", round(x$nalpha), "\n")

## Figure: Hill curve + half-sample rank boxes for the ten leaders.
pdf("manuscript/figures/crime_application.pdf", width = 8.8, height = 3.6)
layout(matrix(1:2, 1), widths = c(1, 1.35))
par(mar = c(4.2, 4.2, 1.4, 0.8))
plot(x$hill_k, x$hill, type = "l", lwd = 1.4, ylim = c(0, .7),
     xlab = "number of upper order statistics k",
     ylab = "Hill estimate")
abline(v = x$nalpha, lty = 2)
par(mar = c(8.6, 4.2, 1.4, 0.6))
boxplot(t(half_mat), names = nm[lead], las = 2, cex.axis = .68,
        outline = FALSE, ylab = "rank in 200 half-samples", ylim = c(0, 40))
abline(h = 20, lty = 3)
dev.off()
cat("WROTE crime outputs\n")
