args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1L) args[1L] else
  "data/real/online_news/OnlineNewsPopularity/OnlineNewsPopularity.csv"
outdir <- if (length(args) >= 2L) args[2L] else "results/real_online_news"
cores <- if (length(args) >= 3L) as.integer(args[3L]) else 1L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

source("code/R/estimator.R")
Rcpp::sourceCpp("code/src/local_hill.cpp")

d <- read.csv(input, check.names = TRUE, stringsAsFactors = FALSE)
names(d) <- trimws(names(d))
stopifnot(nrow(d) == 39644L, ncol(d) == 61L, all(is.finite(d$shares)),
          all(d$shares > 0))

# Continuous, publication-time article descriptors.  The threshold is based on
# covariate support only, never on shares.  Binary/calendar/count variables are
# excluded because the present rank theory assumes continuous marginals.
candidates <- c(
  "n_unique_tokens", "n_non_stop_unique_tokens", "average_token_length",
  "kw_avg_min", "kw_min_avg", "kw_max_avg", "kw_avg_avg",
  paste0("LDA_0", 0:4), "global_subjectivity", "global_sentiment_polarity",
  "global_rate_positive_words", "global_rate_negative_words",
  "avg_positive_polarity", "avg_negative_polarity"
)
stopifnot(all(candidates %in% names(d)))

# Remove image/video-only records and the documented malformed rate records.
valid <- d$n_tokens_content > 0 &
  d$n_unique_tokens >= 0 & d$n_unique_tokens <= 1 &
  d$n_non_stop_unique_tokens >= 0 & d$n_non_stop_unique_tokens <= 1
d <- d[valid, , drop = FALSE]
unique_fraction <- vapply(d[candidates], function(x) length(unique(x)) / length(x), numeric(1))
features <- names(unique_fraction)[unique_fraction >= 0.40]
stopifnot(length(features) >= 10L)
x <- as.matrix(d[features])
y <- d$shares
n <- length(y)

# Tail diagnostics do not choose the screening tuning.
ly <- sort(log(y), decreasing = TRUE)
k_grid <- unique(as.integer(round(exp(seq(log(100), log(5000), length.out = 80)))))
hill <- vapply(k_grid, function(k) mean(ly[seq_len(k)] - ly[k + 1L]), numeric(1))

# Frozen sensitivity grid, including the simulation default.  Smaller scores
# rank first.  Stable ranks across nearby tunings are emphasized over one winner.
grid <- expand.grid(a = c(.20, .25, .30), b = c(.15, .20, .25))
fits <- vector("list", nrow(grid))
for (g in seq_len(nrow(grid))) {
  alpha <- n^(-grid$a[g])
  h <- n^(-grid$b[g]) / 2
  z <- t(vapply(seq_len(ncol(x)), function(j)
    score_coordinate_cpp(x[, j], y, h, alpha, .05), numeric(5)))
  fits[[g]] <- data.frame(a = grid$a[g], b = grid$b[g], alpha = alpha, h = h,
    feature = features, score = z[, "score"], rank = rank(z[, "score"], ties.method = "first"),
    mean_local_n = z[, "mean_local_n"], mean_local_k = z[, "mean_local_k"],
    under_rate = z[, "under_rate"], check.names = FALSE)
}
scores <- do.call(rbind, fits)

# Chronological stability: timedelta is larger for earlier publications.
cut <- median(d$timedelta)
period <- ifelse(d$timedelta > cut, "earlier", "later")
split_fits <- lapply(c("earlier", "later"), function(label) {
  take <- period == label
  nn <- sum(take); alpha <- nn^(-.20); h <- nn^(-.20) / 2
  z <- t(vapply(seq_len(ncol(x)), function(j)
    score_coordinate_cpp(x[take, j], y[take], h, alpha, .05), numeric(5)))
  data.frame(period = label, feature = features, score = z[, "score"],
             rank = rank(z[, "score"], ties.method = "first"))
})
split_scores <- do.call(rbind, split_fits)

# The retained variables are high-cardinality but not perfectly continuous.
# Randomly break exact ties 20 times while preserving every strict ordering.
# This audits sensitivity to the implementation's otherwise stable tie order.
tie_fits <- lapply(seq_len(20L), function(seed) {
  set.seed(83000L + seed)
  xr <- vapply(seq_len(ncol(x)), function(j)
    rank(x[, j], ties.method = "random"), numeric(n))
  alpha <- n^(-.20); h <- n^(-.20) / 2
  z <- t(vapply(seq_len(ncol(xr)), function(j)
    score_coordinate_cpp(xr[, j], y, h, alpha, .05), numeric(5)))
  data.frame(seed = seed, feature = features, score = z[, "score"],
             rank = rank(z[, "score"], ties.method = "first"))
})
tie_scores <- do.call(rbind, tie_fits)

# Sampling uncertainty: 200 seeded random half-samples at the primary rule.
half_one <- function(seed) {
  set.seed(910000L + seed)
  take <- sample.int(n, floor(n / 2), replace = FALSE)
  nn <- length(take); alpha <- nn^(-.20); h <- nn^(-.20) / 2
  z <- t(vapply(seq_len(ncol(x)), function(j)
    score_coordinate_cpp(x[take, j], y[take], h, alpha, .05), numeric(5)))
  data.frame(seed = seed, feature = features, score = z[, "score"],
             rank = rank(z[, "score"], ties.method = "first"))
}
half_list <- if (.Platform$OS.type == "unix" && cores > 1L)
  parallel::mclapply(seq_len(200L), half_one, mc.cores = cores,
                     mc.preschedule = FALSE) else
  lapply(seq_len(200L), half_one)
half_scores <- do.call(rbind, half_list)

rank_summary <- aggregate(rank ~ feature, scores,
  function(v) c(median = median(v), minimum = min(v), maximum = max(v)))
rank_summary <- data.frame(feature = rank_summary$feature,
  median_rank = rank_summary$rank[, "median"], min_rank = rank_summary$rank[, "minimum"],
  max_rank = rank_summary$rank[, "maximum"])
rank_summary <- rank_summary[order(rank_summary$median_rank, rank_summary$max_rank), ]

write.csv(scores, file.path(outdir, "tuning_scores.csv"), row.names = FALSE)
write.csv(split_scores, file.path(outdir, "temporal_stability.csv"), row.names = FALSE)
write.csv(tie_scores, file.path(outdir, "tie_stability.csv"), row.names = FALSE)
write.csv(half_scores, file.path(outdir, "halfsample_stability.csv"), row.names = FALSE)
write.csv(rank_summary, file.path(outdir, "rank_summary.csv"), row.names = FALSE)
write.csv(data.frame(feature = features, unique_fraction = unique_fraction[features]),
          file.path(outdir, "feature_audit.csv"), row.names = FALSE)

pdf(file.path(outdir, "online_news_application.pdf"), width = 10, height = 4.7)
layout(matrix(1:2, nrow = 1), widths = c(1, 1.12))
par(mar = c(4.2, 4.3, 1.2, .7))
plot(k_grid, hill, type = "l", lwd = 2, log = "x", xlab = "Number k of upper shares",
     ylab = "Hill estimate", main = "Upper-tail diagnostic")
abline(v = round(n^(1 - .20)), lty = 2, col = "grey45")
top <- head(rank_summary$feature, 10)
mat <- sapply(top, function(f) scores$rank[scores$feature == f])
display <- c(
  kw_max_avg = "keyword max-average", kw_avg_avg = "keyword average-average",
  LDA_03 = "LDA coordinate 03", global_subjectivity = "global subjectivity",
  LDA_02 = "LDA coordinate 02", avg_positive_polarity = "positive-word polarity",
  global_sentiment_polarity = "global sentiment", LDA_04 = "LDA coordinate 04",
  n_unique_tokens = "lexical diversity", n_non_stop_unique_tokens = "non-stop lexical diversity"
)
colnames(mat) <- unname(display[colnames(mat)])
par(mar = c(4.2, 8.0, 1.2, .7))
boxplot(as.data.frame(mat), horizontal = TRUE, las = 1, xlab = "Rank across 9 tunings",
        main = "Stability of leading descriptors", outline = FALSE, cex.axis = .78)
dev.off()

metadata <- list(
  source_url = "https://archive.ics.uci.edu/static/public/332/online+news+popularity.zip",
  raw_rows = 39644L, analysis_rows = n, excluded_rows = 39644L - n,
  response = "shares", response_summary = summary(y), response_quantiles = quantile(y, c(.9,.95,.99,.995,.999)),
  features = features, tuning_grid = grid, hill_k = k_grid, hill = hill,
  data_sha256 = unname(tools::md5sum(input))
)
saveRDS(list(metadata = metadata, scores = scores, split_scores = split_scores,
             tie_scores = tie_scores, half_scores = half_scores,
             rank_summary = rank_summary),
        file.path(outdir, "online_news_application.rds"))
cat("PASS online-news application", n, "rows", length(features), "features\n")
print(rank_summary)
