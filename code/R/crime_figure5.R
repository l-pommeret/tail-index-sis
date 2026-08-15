## Draft-5 crime figure (manuscript/figures/crime_application.pdf).
## Left/middle panels from results/draft4_real_crime.rds (Hill curve, QQ
## plot); right panel from results/draft5_real_crime_rankfix.rds: the ten
## leaders of the aggregated screen (ranking rule of code/R/rank_rules.R),
## nine per-setting ranks (open circles) and aggregated rank (filled).
## usage: Rscript code/R/crime_figure5.R
x <- readRDS("results/draft4_real_crime.rds")
f <- readRDS("results/draft5_real_crime_rankfix.rds")
stopifnot(identical(x$var_names, f$var_names))
nalpha <- round(x$n * x$alpha)

R <- apply(f$scores, 2L, rank, ties.method = "average")
lead <- order(f$rank_new)[1:10]

pdf("manuscript/figures/crime_application.pdf", width = 8.8, height = 3.6)
layout(matrix(1:3, 1), widths = c(1, .8, 1.25))
par(mar = c(4.2, 4.2, 1.4, 0.8))
plot(x$hill_k, x$hill, type = "l", lwd = 1.4, ylim = c(0, .7),
     xlab = "upper order statistics k", ylab = "Hill estimate")
abline(v = nalpha, lty = 2)
plot(x$qq$x, x$qq$y, pch = 20, cex = .4,
     xlab = "log{(n+1)/i}", ylab = "log of i-th largest response")
fit <- lm(y ~ x, data = x$qq[x$qq$x >= log((x$n + 1) / nalpha), ])
abline(fit, lty = 2)
legend("topleft", bty = "n",
       legend = sprintf("slope %.2f on i <= %d", coef(fit)[2], nalpha))
par(mar = c(8.6, 4.2, 1.4, 0.6))
plot(NA, xlim = c(.5, 10.5), ylim = c(80, 0), xaxt = "n",
     xlab = "", ylab = "rank")
axis(1, at = 1:10, labels = f$var_names[lead], las = 2, cex.axis = .66)
for (i in 1:10)
  points(rep(i, ncol(R)), R[lead[i], ], cex = .7)
points(1:10, f$rank_new[lead], pch = 19)
abline(h = 20, lty = 3)
dev.off()
cat("WROTE manuscript/figures/crime_application.pdf\n")
