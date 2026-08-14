## Draft-2 tables and figures from computed results.
## usage: Rscript code/R/draft2_outputs.R {tuning|sensitivity|comparison}

args <- commandArgs(trailingOnly = TRUE)
what <- args[1L]
dir.create("manuscript/tables", showWarnings = FALSE)
dir.create("manuscript/figures", showWarnings = FALSE)

se <- function(p, n) sqrt(p * (1 - p) / n)

if (what == "tuning") {
  x <- readRDS("results/draft2/tuning.rds")$jobs
  df <- do.call(rbind, lapply(x, function(z)
    data.frame(model = z$model, a = z$a, b = z$b,
               rmax = max(z$active_ranks))))
  ag <- aggregate(df$rmax <= 20, df[c("model", "a", "b")], mean)
  names(ag)[4] <- "sure20"
  write.csv(ag, "results/draft2/tuning_summary.csv", row.names = FALSE)
  ## Common robust region: mean Sure-20 across the four models per cell.
  wide <- reshape(ag, idvar = c("a", "b"), timevar = "model",
                  direction = "wide")
  wide$mean4 <- rowMeans(wide[, grep("sure20", names(wide))])
  wide$min4  <- apply(wide[, grep("sure20\\.", names(wide))], 1, min)
  print(wide[order(-wide$mean4), ][1:8, ])
  ## Heatmap figure: four panels.
  pdf("manuscript/figures/tuning_heatmaps.pdf", width = 8.6, height = 6.6)
  par(mfrow = c(2, 2), mar = c(4.0, 4.0, 2.2, 1.0))
  cols <- hcl.colors(20, "YlOrRd", rev = TRUE)
  aa <- sort(unique(ag$a)); bb <- sort(unique(ag$b))
  for (m in c("M1", "M2", "M3", "M4")) {
    w <- ag[ag$model == m, ]
    zz <- outer(aa, bb, Vectorize(function(x, y)
      w$sure20[w$a == x & w$b == y][1]))
    image(aa, bb, zz, zlim = c(0, 1), col = cols, main = m,
          xlab = "tail exponent a", ylab = "bandwidth exponent b")
    for (i in seq_along(aa)) for (j in seq_along(bb))
      text(aa[i], bb[j], sprintf("%.2f", zz[i, j]), cex = .62)
  }
  dev.off()
  cat("WROTE tuning outputs\n")
} else if (what == "sensitivity") {
  x <- readRDS("results/draft2/sensitivity.rds")$jobs
  df <- do.call(rbind, lapply(x, function(z)
    data.frame(n = z$n, p = z$p, rho = z$rho,
               rmax = max(z$active_ranks))))
  ag <- aggregate(cbind(sure20 = df$rmax <= 20, rmax = df$rmax),
                  df[c("n", "p", "rho")], mean)
  ag$se <- se(ag$sure20, 200)
  write.csv(ag, "results/draft2/sensitivity_summary.csv", row.names = FALSE)
  print(ag)
  ## Sure-d curve at the baseline cell.
  base <- Filter(function(z) z$n == 2000 && z$p == 1000 && z$rho == .25, x)
  dgrid <- c(4, 10, 20, 50, 100, 200, 500)
  sured <- vapply(dgrid, function(d)
    mean(vapply(base, function(z) max(z$active_ranks) <= d, logical(1))),
    numeric(1))
  write.csv(data.frame(d = dgrid, sure = sured),
            "results/draft2/sure_d.csv", row.names = FALSE)
  pdf("manuscript/figures/sure_d.pdf", width = 5.6, height = 4.2)
  par(mar = c(4.2, 4.4, 1.0, 1.0))
  plot(dgrid, sured, type = "b", pch = 19, log = "x", ylim = c(0, 1),
       xlab = "retained dimension d",
       ylab = expression(P(R[max] <= d)))
  abline(v = 4, lty = 3)
  dev.off()
  ## One-factor figure: three panels (n, p, rho).
  pdf("manuscript/figures/sensitivity_factors.pdf", width = 8.6, height = 3.2)
  par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2.0, 0.8))
  w <- ag[ag$p == 1000 & ag$rho == .25, ]; w <- w[order(w$n), ]
  plot(w$n, w$sure20, type = "b", pch = 19, ylim = c(0, 1), log = "x",
       xlab = "n", ylab = "Sure-20", main = "sample size")
  w <- ag[ag$n == 2000 & ag$rho == .25, ]; w <- w[order(w$p), ]
  plot(w$p, w$sure20, type = "b", pch = 19, ylim = c(0, 1), log = "x",
       xlab = "p", ylab = "Sure-20", main = "ambient dimension")
  w <- ag[ag$n == 2000 & ag$p == 1000, ]; w <- w[order(w$rho), ]
  plot(w$rho, w$sure20, type = "b", pch = 19, ylim = c(0, 1),
       xlab = expression(rho), ylab = "Sure-20", main = "dependence")
  dev.off()
  cat("WROTE sensitivity outputs\n")
} else if (what == "comparison") {
  x <- readRDS("results/draft2/comparison.rds")$jobs
  meth <- c(ours = "ranks_ours", yu = "ranks_yu", q95 = "ranks_q95",
            q99 = "ranks_q99")
  rows <- list()
  for (m in c("M1", "M2", "M3", "M4")) {
    w <- Filter(function(z) z$model == m, x)
    for (k in names(meth)) {
      rm_ <- vapply(w, function(z) max(z[[meth[k]]]), numeric(1))
      rows[[paste(m, k)]] <- data.frame(
        model = m, method = k, sure4 = mean(rm_ <= 4),
        sure20 = mean(rm_ <= 20), mean_rmax = mean(rm_))
    }
  }
  out <- do.call(rbind, rows)
  write.csv(out, "results/draft2/comparison_summary.csv", row.names = FALSE)
  print(out, digits = 3)
  ## How often qaSIS(.95) puts inactive 5:8 in its top 4 under M2.
  w <- Filter(function(z) z$model == "M2", x)
  in58 <- mean(vapply(w, function(z)
    sum(z$top20_q95[1:4] %in% 5:8), numeric(1)))
  cat("M2: mean number of U5-U8 in qaSIS(.95) top-4:", in58, "\n")
  ## LaTeX table.
  lab <- c(ours = "Tail-index SIS", yu = "Yoshida--Umezu",
           q95 = "Quantile SIS ($\\tau=.95$)", q99 = "Quantile SIS ($\\tau=.99$)")
  z <- c("\\begin{table}[t]", "\\centering",
    "\\caption{Method comparison at $n=2000$, $p=1000$, $\\rho=0.25$, on 200 common replications per model. Sure-$d$ is the probability that all four active variables rank among the $d$ smallest scores; $\\E(R_{\\max})$ is the mean worst active rank. Monte Carlo standard errors for the probabilities are at most $0.035$.}",
    "\\label{tab:comparison}", "\\begin{tabular}{llrrr}", "\\toprule",
    "Model & Method & Sure-4 & Sure-20 & $\\E(R_{\\max})$\\\\", "\\midrule")
  for (m in c("M1", "M2", "M3", "M4")) {
    for (k in c("ours", "yu", "q95")) {
      w <- out[out$model == m & out$method == k, ]
      z <- c(z, sprintf("%s & %s & %.3f & %.3f & %.1f\\\\",
        if (k == "ours") m else "", lab[k], w$sure4, w$sure20, w$mean_rmax))
    }
    if (m != "M4") z <- c(z, "\\addlinespace")
  }
  z <- c(z, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(z, "manuscript/tables/comparison.tex")
  cat("WROTE comparison outputs\n")
} else stop("unknown")
