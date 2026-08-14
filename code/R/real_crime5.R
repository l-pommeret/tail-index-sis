## Draft-4 real-data application: aggregated tail-index SIS on Communities
## and Crime (Unnormalized).  Median-imputed matrix shared by all methods.
## Aggregation identical to the simulations (code/R/campaign5.R): per-setting
## ranks on the 3x3 block AGRID x BGRID, aggregated by the minimum rank,
## final ranking by increasing aggregated value, ties broken by coordinate.
## usage: Rscript code/R/real_crime5.R OUTPUT.rds

args <- commandArgs(trailingOnly = TRUE)
out_path <- args[1L]
suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/qa_sis.R"); source("code/R/yoshida_umezu.R")

ASTAR <- 0.30; BSTAR <- 0.15
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.10, 0.15, 0.20)

## Load and preprocess exactly as in code/R/real_crime.R (response-blind).
src <- new.env()
lines <- readLines("code/R/real_crime.R")
i1 <- grep("^raw <- read.csv", lines)[1]
i2 <- grep("^score_all", lines)[1] - 1L
eval(parse(text = paste(lines[i1:i2], collapse = "\n")), envir = src)
X <- src$X; y <- src$y; n <- src$n; p <- src$p

score_setting <- function(a, b) {
  al <- n^(-a); hh <- n^(-b) / 2
  vapply(seq_len(p), function(j)
    score_coordinate_cpp(X[, j], y, hh, al, 0.05)[["score"]], numeric(1))
}

S <- matrix(NA_real_, p, 9L); ii <- 0L
setting <- character(9L)
for (a in AGRID) for (b in BGRID) {
  ii <- ii + 1L
  S[, ii] <- score_setting(a, b)
  setting[ii] <- sprintf("a%.2f_b%.2f", a, b)
}
R <- apply(S, 2L, rank, ties.method = "first")
colnames(R) <- setting
amin <- apply(R, 1L, min)
rank_agg <- rank(amin, ties.method = "first")
idx_pub <- which(rep(AGRID, each = 3) == ASTAR & rep(BGRID, 3) == BSTAR)
rank_pub <- R[, idx_pub]

uh <- apply(X, 2L, function(v) rank(v, ties.method = "average") / (n + 1))
yu <- yu_score_matrix(uh, y, k = floor(0.072 * n), h = 1)
rank_yu <- rank(-yu$scores, ties.method = "first")
rank_q <- sapply(c(0.90, 0.95, 0.99), function(tt)
  rank(-qa_sis_scores(uh, y, tau = tt), ties.method = "first"))
colnames(rank_q) <- c("q90", "q95", "q99")

## Tail diagnostics (unchanged from Draft 3).
ys <- sort(y, decreasing = TRUE)
hill_k <- sort(unique(c(seq(20L, 600L, by = 5L), 204L)))
hill <- vapply(hill_k, function(k) mean(log(ys[1:k])) - log(ys[k + 1]),
               numeric(1))
qq <- data.frame(x = log((n + 1) / (1:600)), y = log(ys[1:600]))

saveRDS(list(n = n, p = p, astar = ASTAR, bstar = BSTAR,
             agrid = AGRID, bgrid = BGRID,
             alpha = n^(-ASTAR), h = n^(-BSTAR) / 2,
             var_names = colnames(X),
             setting_ranks = R, rank_agg = rank_agg, rank_pub = rank_pub,
             rank_yu = rank_yu, rank_q = rank_q,
             hill_k = hill_k, hill = hill, qq = qq),
        out_path, compress = "xz")
lead <- order(rank_agg)[1:10]
print(data.frame(variable = colnames(X)[lead], agg = rank_agg[lead],
                 pub = rank_pub[lead], min = apply(R, 1, min)[lead],
                 max9 = apply(R, 1, max)[lead],
                 yu = rank_yu[lead], q95 = rank_q[lead, "q95"]))
cat("WROTE", out_path, "\n")
