## Draft-5 tables and figures.
## usage: Rscript code/R/draft5_outputs.R {tuning|comparison}
## The comparison merges campaign6 (A3, B1, rerun under generate5.R) with
## campaign5 rows for M1 and M3 relabelled A1 and A2: the models, rules and
## seed streams are identical, verified by tests/test_campaign_equivalence.R.
args <- commandArgs(trailingOnly = TRUE)
what <- args[1L]

RULE_MAP <- c("screen seul" = "tail selected",
              "screen 9 minimum" = "tail aggregated",
              "Yoshida-Umezu" = "Yoshida-Umezu",
              "quantile .90" = "quantile .90",
              "quantile .95" = "quantile .95",
              "quantile .99" = "quantile .99")

if (what == "tuning") {
  ag <- read.csv("results/tuning6/summary.csv")
  pdf("manuscript/figures/tuning_heatmaps.pdf", width = 8.6, height = 6.6)
  par(mfrow = c(2, 2), mar = c(4.0, 4.0, 2.2, 1.0))
  cols <- hcl.colors(20, "YlOrRd", rev = TRUE)
  aa <- sort(unique(ag$a)); bb <- sort(unique(ag$b))
  for (m in c("A1", "A2", "A3", "B1")) {
    w <- ag[ag$model == m, ]
    zz <- outer(aa, bb, Vectorize(function(x, y)
      w$sure20[w$a == x & w$b == y][1]))
    image(aa, bb, zz, zlim = c(0, 1), col = cols, main = m,
          xlab = "tail exponent a", ylab = "bandwidth exponent b")
    for (i in seq_along(aa)) for (j in seq_along(bb))
      text(aa[i], bb[j], sprintf("%.2f", zz[i, j]), cex = .6)
    rect(0.275, 0.075, 0.425, 0.225, border = "grey30", lwd = 1.4, lty = 2)
  }
  dev.off()
  w <- reshape(ag, idvar = c("a", "b"), timevar = "model",
               direction = "wide")
  w$mean4 <- rowMeans(w[, grep("sure20", names(w))])
  print(w[order(-w$mean4), ][1:10, ], digits = 3)
  cat("WROTE tuning heatmaps\n")
} else if (what == "comparison") {
  c6 <- read.csv("results/campaign6/summary.csv")
  c5 <- read.csv("results/campaign5/summary.csv")
  c5 <- c5[c5$model %in% c("M1", "M3"), ]
  c5$model <- ifelse(c5$model == "M1", "A1", "A2")
  c5$rule <- RULE_MAP[c5$rule]
  ## align columns (campaign5 has identical metric columns)
  common <- intersect(names(c6), names(c5))
  all <- rbind(c6[, common], c5[, common])
  write.csv(all, "results/draft5_comparison_summary.csv", row.names = FALSE)
  lab <- c("tail selected" = "Tail-index SIS, selected $(a^\\star,b^\\star)$",
           "tail aggregated" = "Tail-index SIS, aggregated",
           "Yoshida-Umezu" = "Yoshida--Umezu",
           "quantile .90" = "Quantile SIS ($\\tau=.90$)",
           "quantile .95" = "Quantile SIS ($\\tau=.95$)",
           "quantile .99" = "Quantile SIS ($\\tau=.99$)")
  for (pp in c(500, 1000, 2000)) {
    z <- c("\\begin{table}[t]", "\\centering",
      sprintf("\\caption{Method comparison at $n=2000$, $p=%d$, $\\rho=0.25$, on the same 1{,}000 replications per model. Sure-$d$ is the probability that all four tail-index-active variables rank among the $d$ smallest scores ($d=4$ is exact recovery); $\\mathbb{E}(R_{\\max})$ and $\\mathrm{Med}(R_{\\max})$ are the mean and median of the worst active rank $R_{\\max}=\\max_{j\\in\\A_\\gamma}\\mathrm{rank}(j)$. Standard errors of the probabilities are at most $0.016$.}", pp),
      sprintf("\\label{tab:cmp%d}", pp), "\\begin{tabular}{llrrrr}",
      "\\toprule",
      "Model & Method & Sure-4 & Sure-20 & $\\mathbb{E}(R_{\\max})$ & $\\mathrm{Med}(R_{\\max})$\\\\",
      "\\midrule")
    for (m in c("A1", "A2", "A3", "B1")) {
      for (k in names(lab)) {
        w <- all[all$p == pp & all$model == m & all$rule == k, ]
        z <- c(z, sprintf("%s & %s & %.3f & %.3f & %.1f & %.0f\\\\",
          ifelse(k == "tail selected", m, ""), lab[k],
          w$sure4, w$sure20, w$ermax, w$medmax))
      }
      if (m != "B1") z <- c(z, "\\addlinespace")
    }
    z <- c(z, "\\bottomrule", "\\end{tabular}", "\\end{table}")
    writeLines(z, sprintf("manuscript/tables/cmp%d.tex", pp))
  }
  print(all[all$model %in% c("A3", "B1") &
            all$rule %in% c("tail selected", "tail aggregated",
                            "quantile .95"), ], digits = 3)
  cat("WROTE 3 comparison tables\n")
} else stop("unknown")
