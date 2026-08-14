## Draft-3 tables and figures.
## usage: Rscript code/R/draft3_outputs.R {tuning|comparison|real}

args <- commandArgs(trailingOnly = TRUE)
what <- args[1L]

if (what == "tuning") {
  ## Merge draft-2 heatmaps for M1/M3/M4 with the draft-3 M2 heatmap.
  files <- list.files("results/draft2/tuning_cells", full.names = TRUE)
  files <- files[!grepl("_M2_", files)]
  rows <- do.call(rbind, lapply(files, function(f) {
    x <- readRDS(f)
    do.call(rbind, lapply(x, function(z)
      data.frame(model = z$model, a = z$a, b = z$b,
                 rmax = max(z$active_ranks))))
  }))
  f2 <- list.files("results/draft3/tuning_m2_cells", full.names = TRUE)
  rows2 <- do.call(rbind, lapply(f2, function(f) {
    x <- readRDS(f)
    do.call(rbind, lapply(x, function(z)
      data.frame(model = "M2", a = z$a, b = z$b,
                 rmax = max(z$active_ranks))))
  }))
  ag <- aggregate(rmax <= 20 ~ model + a + b, rbind(rows, rows2), mean)
  names(ag)[4] <- "sure20"
  write.csv(ag, "results/draft3/tuning_summary.csv", row.names = FALSE)
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
      text(aa[i], bb[j], sprintf("%.2f", zz[i, j]), cex = .6)
  }
  dev.off()
  w <- ag[ag$a == .30 & ag$b %in% c(.05, .10, .15, .20), ]
  print(reshape(w, idvar = "b", timevar = "model", direction = "wide"))
  cat("WROTE tuning outputs\n")
} else if (what == "comparison") {
  files <- list.files("results/draft3/comparison_cells", full.names = TRUE)
  rows <- list(); taus <- c("q900", "q950", "q975", "q990")
  comp <- list()
  for (f in files) {
    x <- readRDS(f); m <- x[[1]]$model
    meth <- list(ours = lapply(x, `[[`, "ranks_ours"),
                 yu = lapply(x, `[[`, "ranks_yu"))
    for (tt in taus) meth[[tt]] <- lapply(x, function(z) z$q[[tt]]$ranks)
    for (k in names(meth)) {
      rm_ <- sapply(meth[[k]], max)
      rows[[paste(m, k)]] <- data.frame(model = m, method = k,
        sure4 = mean(rm_ <= 4), sure20 = mean(rm_ <= 20),
        ermax = round(mean(rm_), 2))
    }
    ## top-4 composition (A_gamma = 1:4, A_scale = 5:8) for every method
    t4 <- list(ours = lapply(x, `[[`, "top4_ours"),
               yu = lapply(x, `[[`, "top4_yu"))
    for (tt in taus) t4[[tt]] <- lapply(x, function(z) z$q[[tt]]$top4)
    for (k in names(t4)) {
      comp[[paste(m, k)]] <- data.frame(model = m, method = k,
        top4_gamma = mean(sapply(t4[[k]], function(v) sum(v %in% 1:4))),
        top4_scale = mean(sapply(t4[[k]], function(v) sum(v %in% 5:8))))
    }
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  cc <- do.call(rbind, comp); rownames(cc) <- NULL
  write.csv(out, "results/draft3/comparison_summary.csv", row.names = FALSE)
  write.csv(cc, "results/draft3/top4_composition.csv", row.names = FALSE)
  print(out, digits = 3); print(cc, digits = 3)
  ## tau-sensitivity figure
  pdf("manuscript/figures/tau_sensitivity.pdf", width = 5.8, height = 4.4)
  par(mar = c(4.2, 4.4, 1.0, 1.0))
  tv <- c(.90, .95, .975, .99)
  plot(range(tv), c(0, 1), type = "n", xlab = expression(tau),
       ylab = "Sure-20 of quantile SIS")
  pchs <- c(M1 = 19, M2 = 17, M3 = 15, M4 = 18)
  for (m in c("M1", "M2", "M3", "M4")) {
    v <- sapply(taus, function(tt) out$sure20[out$model == m & out$method == tt])
    lines(tv, v, type = "b", pch = pchs[m])
    ours <- out$sure20[out$model == m & out$method == "ours"]
    points(0.9885, ours, pch = pchs[m], col = "grey55", xpd = NA)
  }
  legend("bottomleft", legend = c("M1", "M2", "M3", "M4"),
         pch = pchs, lty = 1, bty = "n")
  dev.off()
  cat("WROTE comparison outputs\n")
} else if (what == "real") {
  x <- readRDS("results/draft3/real_crime.rds")
  p <- x$p; nm <- x$var_names
  primary_rank <- rank(x$primary_scores, ties.method = "first")
  lead <- order(x$primary_scores)[1:10]
  grid_mat <- vapply(x$grid_ranks, function(r) r[lead], numeric(10))
  half_mat <- vapply(x$half_ranks, function(r) r[lead], numeric(10))
  tie_mat <- vapply(x$tie_ranks, function(r) r[lead], numeric(10))
  q1 <- function(v, pr) as.numeric(quantile(v, pr, type = 1))
  summ <- data.frame(variable = nm[lead], primary = primary_rank[lead],
    grid_med = apply(grid_mat, 1, q1, .5), grid_max = apply(grid_mat, 1, max),
    half_med = apply(half_mat, 1, q1, .5),
    half_q10 = apply(half_mat, 1, q1, .1),
    half_q90 = apply(half_mat, 1, q1, .9),
    top5_freq = rowMeans(half_mat <= 5), tie_max = apply(tie_mat, 1, max),
    yu_rank = rank(-x$yu_scores, ties.method = "first")[lead],
    q95_rank = rank(-x$q95_scores, ties.method = "first")[lead])
  write.csv(summ, "results/draft3/real_crime_leaders.csv", row.names = FALSE)
  print(summ, digits = 3)
  pdf("manuscript/figures/crime_application.pdf", width = 8.8, height = 3.6)
  layout(matrix(1:3, 1), widths = c(1, .8, 1.25))
  par(mar = c(4.2, 4.2, 1.4, 0.8))
  plot(x$hill_k, x$hill, type = "l", lwd = 1.4, ylim = c(0, .7),
       xlab = "upper order statistics k", ylab = "Hill estimate")
  abline(v = x$nalpha, lty = 2)
  plot(x$qq$x, x$qq$y, pch = 20, cex = .4,
       xlab = "log{(n+1)/i}", ylab = "log of i-th largest response")
  ## slope fitted over the working range i <= n*alpha (~204): estimates gamma
  fit <- lm(y ~ x, data = x$qq[x$qq$x >= log((x$n + 1) / round(x$nalpha)), ])
  abline(fit, lty = 2)
  legend("topleft", bty = "n",
         legend = sprintf("slope %.2f on i <= %d", coef(fit)[2],
                          round(x$nalpha)))
  par(mar = c(8.6, 4.2, 1.4, 0.6))
  boxplot(t(half_mat), names = nm[lead], las = 2, cex.axis = .66,
          outline = FALSE, ylab = "rank in 200 half-samples", ylim = c(0, 45))
  abline(h = 20, lty = 3)
  dev.off()
  cat("WROTE real outputs; n =", x$n, "p =", p,
      "2nalphah =", round(2 * x$n * x$alpha * x$h), "\n")
} else stop("unknown")
