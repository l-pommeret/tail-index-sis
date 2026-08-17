## Draft-6 tables and figures (article_light_2).
## usage: Rscript code/R/draft6_outputs.R {comparison|tuning|yu}
##
## Unlike draft5_outputs.R, the comparison reads a single source,
## results/campaign8/summary.csv: campaign8 recomputes every rule on every
## model with the permutation-invariant ranking rules, so no merging with
## campaign5/6/7 is required and no relabelling is involved.
args <- commandArgs(trailingOnly = TRUE)
what <- if (length(args)) args[1L] else "comparison"

## Target article directory: second argument, or $D6_ARTICLE, defaulting to
## article_light_2.  Both live versions share the same tables and figure.
ART <- if (length(args) >= 2) args[2L] else
  Sys.getenv("D6_ARTICLE", "article_light_2")
OUT_TAB <- file.path(ART, "tables")
OUT_FIG <- file.path(ART, "figures")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

if (what == "tuning") {
  ag <- read.csv("results/tuning6/summary.csv")
  pdf(file.path(OUT_FIG, "tuning_heatmaps.pdf"), width = 8.6, height = 6.6)
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
    ## dashed: the 3x3 aggregation block N9; solid: the selected cell.
    ## Drawn as borders, not markers, so the printed Sure-20 stays legible.
    rect(0.275, 0.075, 0.425, 0.225, border = "grey30", lwd = 1.4, lty = 2)
    rect(0.325, 0.125, 0.375, 0.175, border = "black", lwd = 2.6)
  }
  dev.off()
  cat("WROTE tuning heatmaps (selected cell marked at a=0.35, b=0.15)\n")

} else if (what == "yu") {
  tab <- read.csv("results/yu_tuning8/yu_tuning_grid8.csv")
  sel <- read.csv("results/yu_tuning8/yu_selected.csv")
  z <- c("\\begin{table}[t]", "\\centering",
    "\\caption{Tuning grid for the conditional-Pickands screen of \\citet{yoshidaUmezu2026}, Sure-20 over $200$ common-random-number replications per model at $n=2000$, $p=1000$, $\\rho=0.25$ --- the same datasets used to tune the tail-index screen. The retained cell for each model is shown in bold. The grid brackets the bandwidth range and the intermediate-sequence range examined in their own sensitivity study.}",
    "\\label{tab:yutuning}", "\\begin{tabular}{lrrrrr}", "\\toprule",
    "Model & $h$ & $k$ & $k/n$ & Sure-4 & Sure-20\\\\", "\\midrule")
  for (m in c("A1", "A2", "A3", "B1")) {
    w <- tab[tab$model == m, ]
    w <- w[order(w$h, w$k), ]
    bh <- sel$h[sel$model == m]; bk <- sel$k[sel$model == m]
    for (i in seq_len(nrow(w))) {
      star <- w$h[i] == bh && w$k[i] == bk
      fmt <- if (star) "\\textbf{%.3f}" else "%.3f"
      z <- c(z, sprintf(paste0("%s & %.1f & %d & %.3f & ", fmt, " & ", fmt,
                               "\\\\"),
        ifelse(i == 1L, m, ""), w$h[i], w$k[i], w$k_fraction[i],
        w$sure4[i], w$sure20[i]))
    }
    if (m != "B1") z <- c(z, "\\addlinespace")
  }
  z <- c(z, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(z, file.path(OUT_TAB, "yutuning.tex"))
  cat("WROTE YU tuning table\n"); print(sel, row.names = FALSE)

} else if (what == "comparison") {
  all <- read.csv("results/campaign8/summary.csv")
  write.csv(all, "results/draft6_comparison_summary.csv", row.names = FALSE)
  lab <- c(
    "tail selected"       = "Tail-index SIS, selected $(a^\\star,b^\\star)$",
    "tail aggregated"     = "Tail-index SIS, aggregated",
    "Yoshida-Umezu tuned" = "Yoshida--Umezu, tuned per model",
    "Yoshida-Umezu paper" = "Yoshida--Umezu, paper tuning",
    "quantile .90"        = "Quantile SIS ($\\tau=.90$)",
    "quantile .95"        = "Quantile SIS ($\\tau=.95$)",
    "quantile .99"        = "Quantile SIS ($\\tau=.99$)")
  for (pp in c(500, 1000, 2000)) {
    z <- c("\\begin{table}[t]", "\\centering",
      sprintf("\\caption{Method comparison at $n=2000$, $p=%d$, $\\rho=0.25$, on the same 1{,}000 replications per model. Sure-$d$ is the probability that all four tail-index-active variables rank among the $d$ smallest scores ($d=4$ is exact recovery); $\\mathbb{E}(R_{\\max})$ and $\\mathrm{Med}(R_{\\max})$ are the mean and median of the worst active rank $R_{\\max}=\\max_{j\\in\\A_\\gamma}\\mathrm{rank}(j)$. The tail-index screen is reported at one $(a^\\star,b^\\star)$ common to all four models; the conditional-Pickands screen is reported both at its paper tuning and at the cell of its own grid retained separately for each model (Table~\\ref{tab:yutuning}). Standard errors of the probabilities are at most $0.016$.}", pp),
      sprintf("\\label{tab:cmp%d}", pp), "\\begin{tabular}{llrrrr}",
      "\\toprule",
      "Model & Method & Sure-4 & Sure-20 & $\\mathbb{E}(R_{\\max})$ & $\\mathrm{Med}(R_{\\max})$\\\\",
      "\\midrule")
    for (m in c("A1", "A2", "A3", "B1")) {
      for (k in names(lab)) {
        w <- all[all$p == pp & all$model == m & all$rule == k, ]
        if (nrow(w) != 1L) stop("missing cell: ", m, " p=", pp, " ", k)
        z <- c(z, sprintf("%s & %s & %.3f & %.3f & %.1f & %.0f\\\\",
          ifelse(k == "tail selected", m, ""), lab[k],
          w$sure4, w$sure20, w$ermax, w$medmax))
      }
      if (m != "B1") z <- c(z, "\\addlinespace")
    }
    z <- c(z, "\\bottomrule", "\\end{tabular}", "\\end{table}")
    writeLines(z, file.path(OUT_TAB, sprintf("cmp%d.tex", pp)))
  }
  ## Table: composition of the leading positions on B1.  Regenerated here
  ## as well; it was previously carried over by hand and so kept the
  ## numbers of an earlier campaign.
  z <- c("\\begin{table}[t]", "\\centering",
    "\\caption{Model B1: average number of the twenty scale proxies among the top 4 and top 24 positions of each ranking, by dimension (same $1{,}000$ replications as Tables~\\ref{tab:cmp500}--\\ref{tab:cmp2000}). A perfectly specific screen scores $0$; ranking the proxies first scores $4$ and $20$.}",
    "\\label{tab:b1comp}", "\\begin{tabular}{lrrrrrr}", "\\toprule",
    " & \\multicolumn{2}{c}{$p=500$} & \\multicolumn{2}{c}{$p=1000$}",
    " & \\multicolumn{2}{c}{$p=2000$}\\\\",
    "Method & top 4 & top 24 & top 4 & top 24 & top 4 & top 24\\\\",
    "\\midrule")
  for (k in names(lab)) {
    w <- lapply(c(500, 1000, 2000), function(pp)
      all[all$p == pp & all$model == "B1" & all$rule == k, ])
    if (any(vapply(w, nrow, 0L) != 1L)) stop("missing B1 cell: ", k)
    z <- c(z, sprintf("%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f\\\\",
      lab[k],
      w[[1]]$top4_scale, w[[1]]$top24_scale,
      w[[2]]$top4_scale, w[[2]]$top24_scale,
      w[[3]]$top4_scale, w[[3]]$top24_scale))
  }
  z <- c(z, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(z, file.path(OUT_TAB, "b1comp.tex"))

  ## What the competitor's grid search bought it.
  d <- merge(all[all$rule == "Yoshida-Umezu paper",
                 c("model", "p", "sure4", "sure20", "ermax")],
             all[all$rule == "Yoshida-Umezu tuned",
                 c("model", "p", "sure4", "sure20", "ermax")],
             by = c("model", "p"), suffixes = c("_paper", "_tuned"))
  cat("\nEffect of the grid search on Yoshida--Umezu:\n")
  print(d[order(d$p, d$model), ], row.names = FALSE, digits = 3)
  cat("\nWROTE 3 comparison tables\n")
} else stop("unknown target: ", what)
