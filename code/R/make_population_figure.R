args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else "manuscript/figures/population_profiles.pdf"
u <- seq(0, 1, length.out = 401)
pdf(out, width = 9.2, height = 3.8)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 1.5, .8))
plot(u, .5 * exp(-u), type = "l", lwd = 2, ylim = c(.16, .53),
     xlab = expression(u), ylab = expression(xi[j](u)),
     main = "Detectable monotone activity")
abline(h = .5, lwd = 2, lty = 2)
legend("bottomleft", c("active coordinate", "inactive coordinate"),
       lty = c(1, 2), lwd = 2, bty = "n")

plot(u, rep(1.2, length(u)), type = "l", lwd = 2, ylim = c(.84, 1.23),
     xlab = expression(u), ylab = expression(xi[j](u)),
     main = "Active but envelope-invisible")
text(.5, 1.145, "flat projected profile", cex = .9)
text(.5, .99, expression(gamma(u[1],u[2]) == 1.2-(u[1]-u[2])^2), cex = .9)
text(.5, .91, "every fibre meets the diagonal maximizer", cex = .82)
dev.off()
cat("WROTE", out, "\n")
