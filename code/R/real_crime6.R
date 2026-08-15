## Crime application under the permutation-invariant ranking rule
## (audit of the column-index tie-break; see code/R/rank_rules.R).
## Recomputes the aggregated and selected rankings from the same score
## matrix construction as real_crime5.R, reports both OLD (index tie-break)
## and NEW (invariant) leader tables, and their differences.
## Dedicated tie-break seed: 973000019 (real data, no replication index).
## usage: Rscript code/R/real_crime6.R OUTPUT.rds
args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1) args[1] else
  "results/draft5_real_crime_rankfix.rds"
suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/rank_rules.R")

ASTAR <- 0.30; BSTAR <- 0.15
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.10, 0.15, 0.20)

src <- new.env()
lines <- readLines("code/R/real_crime.R")
i1 <- grep("^raw <- read.csv", lines)[1]
i2 <- grep("^score_all", lines)[1] - 1L
eval(parse(text = paste(lines[i1:i2], collapse = "\n")), envir = src)
X <- src$X; y <- src$y; n <- src$n; p <- src$p

S <- matrix(NA_real_, p, 9L); ii <- 0L
for (a in AGRID) for (b in BGRID) {
  ii <- ii + 1L
  al <- n^(-a); hh <- n^(-b) / 2
  S[, ii] <- vapply(seq_len(p), function(j)
    score_coordinate_cpp(X[, j], y, hh, al, 0.05)[["score"]], numeric(1))
}
idx_pub <- which(rep(AGRID, each = 3) == ASTAR & rep(BGRID, 3) == BSTAR)

u <- tiebreak_u(p, 973000019L)
o_old <- order_agg_old(S);       r_old <- integer(p); r_old[o_old] <- 1:p
o_new <- order_agg_new(S, u);    r_new <- integer(p); r_new[o_new] <- 1:p
key <- agg_key(setting_ranks_avg(S))
ts <- tie_stats(S, integer(0))

lead <- function(o, r, k = 10L) data.frame(
  variable = colnames(X)[o[1:k]], rank = r[o[1:k]],
  min9 = key[o[1:k], "min"], med9 = key[o[1:k], "med"],
  max9 = key[o[1:k], "max"])
cat("== OLD (index tie-break) leaders ==\n"); print(lead(o_old, r_old))
cat("== NEW (invariant) leaders ==\n");       print(lead(o_new, r_new))
moved <- which(r_old != r_new)
cat(sprintf("\n%d of %d variables change rank; max |shift| = %d\n",
            length(moved), p,
            if (length(moved)) max(abs(r_old - r_new)) else 0L))
if (length(moved))
  print(data.frame(variable = colnames(X)[moved], old = r_old[moved],
                   new = r_new[moved])[order(pmin(r_old[moved],
                                                  r_new[moved])), ][1:min(15,
                                                  length(moved)), ])
cat(sprintf("amin ties: %d of %d coordinates share their aggregated value\n",
            ts[["tied_coords_amin"]], p))
saveRDS(list(n = n, p = p, var_names = colnames(X), scores = S,
             key = key, rank_old = r_old, rank_new = r_new,
             tiebreak_seed = 973000019L, tie_stats = ts),
        out_path, compress = "xz")
cat("WROTE", out_path, "\n")
