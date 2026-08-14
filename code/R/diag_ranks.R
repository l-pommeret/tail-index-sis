## Why does a Sure-d curve saturate fast or slowly?  Sure-d is the c.d.f. of
## Rmax = max_j rank(active j), so its shape is the shape of that distribution.
## This prints, per cell and method: the quantiles of Rmax, the per-coordinate
## mean rank (is one active coordinate structurally weaker?), and the share of
## replications that fail catastrophically.
## usage: Rscript code/R/diag_ranks.R [CELLDIR] [pattern]

args <- commandArgs(trailingOnly = TRUE)
celldir <- if (length(args) >= 1L) args[1L] else "results/grid/comparison_cells"
pat <- if (length(args) >= 2L) args[2L] else "\\.rds$"
METHODS <- c("ours", "yu", "q900", "q950", "q975", "q990")
ranks_of <- function(z, k)
  if (k == "ours") z$ranks_ours else if (k == "yu") z$ranks_yu else z$q[[k]]$ranks

for (f in sort(list.files(celldir, pattern = pat, full.names = TRUE))) {
  x <- readRDS(f); z1 <- x[[1]]
  cat(sprintf("\n=== %s  (n=%d p=%d rho=%.2f, %d reps)\n", z1$model, z1$n,
              z1$p, z1$rho, length(x)))
  cat(sprintf("%-22s %5s %5s %5s %5s %5s %5s | %6s %6s %6s %6s | %5s %5s\n",
              "methode", "q10", "q25", "med", "q75", "q90", "max",
              "r(X1)", "r(X2)", "r(X3)", "r(X4)", ">20", ">100"))
  for (k in METHODS) {
    R <- vapply(x, function(z) as.numeric(ranks_of(z, k)), numeric(4))  # 4 x reps
    rmax <- apply(R, 2L, max)
    qq <- quantile(rmax, c(.1, .25, .5, .75, .9, 1), type = 1)
    cat(sprintf("%-22s %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f | %6.1f %6.1f %6.1f %6.1f | %5.2f %5.2f\n",
                k, qq[1], qq[2], qq[3], qq[4], qq[5], qq[6],
                mean(R[1, ]), mean(R[2, ]), mean(R[3, ]), mean(R[4, ]),
                mean(rmax > 20), mean(rmax > 100)))
  }
  ## which coordinate carries the worst rank, for the proposed screen
  R <- vapply(x, function(z) as.numeric(z$ranks_ours), numeric(4))
  who <- apply(R, 2L, which.max)
  cat("  argmax du pire rang (Tail-index SIS) :",
      paste(sprintf("X%d=%.0f%%", 1:4, 100 * tabulate(who, 4) / length(who)),
            collapse = "  "), "\n")
}
