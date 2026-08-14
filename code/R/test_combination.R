## Combining the tail-index screen with the quantile screen, evaluated on the
## grid cells already computed -- no new simulation, since every cell stores the
## rank of each active coordinate under every method.
##
##   min-rank  a coordinate is kept as soon as EITHER screen ranks it in its
##             top d: the union rule, recall-oriented.  It keeps at most 2d
##             coordinates, so the size-fair comparison is against a single
##             screen at 2d, reported alongside.
##   max-rank  a coordinate is kept only if BOTH screens rank it in their top d:
##             the intersection rule, specificity-oriented.  A scale-only
##             variable is ranked high by the quantile screen and low by the
##             tail-index screen, so the intersection is precisely what a
##             scale-only variable fails.
##
## usage: Rscript code/R/test_combination.R [CELLDIR] [TAU_KEY]

args <- commandArgs(trailingOnly = TRUE)
celldir <- if (length(args) >= 1L) args[1L] else "results/grid/comparison_cells"
tkey <- if (length(args) >= 2L) args[2L] else "q950"
DS <- c(4, 10, 20, 30, 50, 100)

files <- sort(list.files(celldir, pattern = "\\.rds$", full.names = TRUE))
rows <- list()
for (f in files) {
  x <- readRDS(f); z1 <- x[[1]]
  rt <- vapply(x, function(z) as.numeric(z$ranks_ours), numeric(4))
  rq <- vapply(x, function(z) as.numeric(z$q[[tkey]]$ranks), numeric(4))
  rmax <- list(
    ours   = apply(rt, 2L, max),
    qsis   = apply(rq, 2L, max),
    union  = apply(pmin(rt, rq), 2L, max),
    inter  = apply(pmax(rt, rq), 2L, max))
  for (k in names(rmax)) {
    r <- data.frame(model = z1$model, n = z1$n, p = z1$p, rho = z1$rho,
                    rule = k, reps = length(x))
    for (d in DS) r[[paste0("sure", d)]] <- mean(rmax[[k]] <= d)
    ## size-fair column: a union kept at d holds up to 2d coordinates
    r$sure20_fair <- mean(rmax[[k]] <= if (k == "union") 10 else 20)
    rows[[length(rows) + 1L]] <- r
  }
}
S <- do.call(rbind, rows)
write.csv(S, file.path(dirname(celldir), "combination_summary.csv"),
          row.names = FALSE)

fm <- function(v) formatC(v, format = "f", digits = 3)
cat(sprintf("\nRegles combinees, quantile SIS a %s, %d cellules\n\n", tkey,
            length(files)))
cat(sprintf("%-6s %-7s %s %12s\n", "modele", "regle",
            paste(sprintf("%8s", paste0("S", DS)), collapse = ""),
            "S20 equitable"))
for (m in c("M1", "M2", "M3", "M4")) {
  w <- S[S$model == m, ]
  if (!nrow(w)) next
  for (k in c("ours", "qsis", "union", "inter")) {
    v <- vapply(paste0("sure", DS), function(cc) mean(w[[cc]][w$rule == k]),
                numeric(1))
    cat(sprintf("%-6s %-7s %s %12s\n", m, k,
                paste(sprintf("%8s", fm(v)), collapse = ""),
                fm(mean(w$sure20_fair[w$rule == k]))))
  }
  cat("\n")
}
