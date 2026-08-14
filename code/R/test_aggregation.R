## Two questions at once, sharing the expensive part (the score passes):
##   (a) does a finer tuning grid beat a coarser one, and does a wider grid beat
##       a narrower one?  Denser grids add correlated rankings, wider grids add
##       diverse but individually worse ones.
##   (b) is the median the right way to aggregate ranks across tunings?
##
## Grids: coarse 3x3 over the good region, fine 5x5 over the SAME region,
## wide 5x5 over the paper's full tuning range (which includes settings measured
## to be poor: SNR 1.36 at b=0.20 against 3.10 at b=0.10).
## Aggregators: median, mean, 20% trimmed mean, minimum, maximum, geometric mean.
## Each is evaluated standalone and behind the quantile-screen prefilter.
##
## usage: KAPPA=0.20 Rscript code/R/test_aggregation.R CORES CELL.rds [...]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R"); source("code/R/qa_sis.R")

args <- commandArgs(trailingOnly = TRUE)
CORES <- as.integer(args[1L]); cells <- args[-1L]
dir.create("results/aggregation", showWarnings = FALSE, recursive = TRUE)
EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
D1 <- 25L

GRIDS <- list(
  coarse = list(a = c(0.30, 0.35, 0.40), b = c(0.05, 0.10, 0.15)),
  fine   = list(a = seq(0.300, 0.400, by = 0.025),
                b = seq(0.050, 0.150, by = 0.025)),
  wide   = list(a = seq(0.25, 0.45, by = 0.05),
                b = seq(0.00, 0.20, by = 0.05)))
AGGS <- list(
  mediane = function(R) apply(R, 1L, median),
  moyenne = function(R) rowMeans(R),
  tronq20 = function(R) apply(R, 1L, mean, trim = 0.2),
  minimum = function(R) apply(R, 1L, min),
  maximum = function(R) apply(R, 1L, max),
  geom    = function(R) exp(rowMeans(log(R))))

score_vec <- function(z, y, a, b) {
  n <- nrow(z); p <- ncol(z); h <- n^(-b) / 2; alpha <- n^(-a)
  vapply(seq_len(p), function(j)
    score_coordinate_cpp(z[, j], y, h, alpha, EPS), numeric(5))["score", ]
}

one <- function(z0) {
  n <- z0$n; p <- z0$p
  d <- simulate_dataset3(n, p, z0$rho, z0$model, z0$seed)
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (n + 1))
  ord_q <- order(-qa_sis_scores(uh, d$y, tau = 0.95), seq_len(p))
  surv <- ord_q[seq_len(D1)]
  base <- score_vec(d$z, d$y, ASTAR, BSTAR)
  out <- list(seul = max(match(1:4, order(base, seq_len(p), na.last = TRUE))),
              quantile = max(match(1:4, ord_q)))
  ## cache every tuning used by any grid, so no score is computed twice
  need <- unique(do.call(rbind, lapply(GRIDS, function(g)
    expand.grid(a = g$a, b = g$b))))
  key <- sprintf("%.3f_%.3f", need$a, need$b)
  ranks <- list()
  for (i in seq_len(nrow(need))) {
    s <- score_vec(d$z, d$y, need$a[i], need$b[i])
    ranks[[key[i]]] <- rank(s, ties.method = "first", na.last = TRUE)
  }
  for (gname in names(GRIDS)) {
    g <- GRIDS[[gname]]
    gg <- expand.grid(a = g$a, b = g$b)
    R <- do.call(cbind, ranks[sprintf("%.3f_%.3f", gg$a, gg$b)])
    for (aname in names(AGGS)) {
      agg <- AGGS[[aname]](R)
      o <- order(agg, seq_len(p))
      out[[paste(gname, aname, sep = ".")]] <-
        max(match(1:4, o))
      out[[paste(gname, aname, "pipe", sep = ".")]] <-
        if (!all(1:4 %in% surv)) Inf
        else max(match(1:4, surv[order(agg[surv], surv)]))
    }
  }
  c(list(model = z0$model, n = n, p = p, rho = z0$rho), out)
}

res <- list()
for (f in cells) {
  x <- readRDS(f)
  o <- parallel::mclapply(x, one, mc.cores = CORES, mc.preschedule = FALSE)
  bad <- vapply(o, function(z) inherits(z, "try-error") || is.null(z), logical(1))
  if (any(bad)) stop(f, ": ", sum(bad), " failed: ", as.character(o[bad][[1]]))
  res[[f]] <- o
  cat(format(Sys.time(), "%H:%M:%S"), basename(f), "done\n"); flush.console()
}
saveRDS(res, file.path("results/aggregation", sprintf("raw_n%d_p%d.rds", res[[1]][[1]]$n, res[[1]][[1]]$p)), compress = "xz")

fm <- function(v) formatC(v, format = "f", digits = 3)
nm <- setdiff(names(res[[1]][[1]]), c("model", "n", "p", "rho"))
cat(sprintf("\nGrilles de reglages et agregateurs (%d, %d et %d points)\n\n",
            length(GRIDS$coarse$a) * length(GRIDS$coarse$b),
            length(GRIDS$fine$a) * length(GRIDS$fine$b),
            length(GRIDS$wide$a) * length(GRIDS$wide$b)))
for (f in cells) {
  o <- res[[f]]; z1 <- o[[1]]
  cat(sprintf("=== %s  n=%d p=%d rho=%.2f\n", z1$model, z1$n, z1$p, z1$rho))
  cat(sprintf("%-26s %8s %8s %8s\n", "regle", "Sure-4", "Sure-10", "Sure-20"))
  for (k in nm) {
    v <- vapply(o, function(z) z[[k]], numeric(1))
    cat(sprintf("%-26s %8s %8s %8s\n", k, fm(mean(v <= 4)), fm(mean(v <= 10)),
                fm(mean(v <= 20))))
  }
  cat("\n")
}
