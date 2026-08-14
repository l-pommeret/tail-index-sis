#!/usr/bin/env Rscript

# Reconciliation benchmark for code_alex/test_grande_dimension.R.
#
# The Alex estimator below is reimplemented from the following line ranges of
# that script, without sourcing it (which would launch its 100-replication run):
#   stable marginal k: lines 416--532;
#   exceedances:       lines 535--573;
#   kernel score:      lines 72--248.
# The article estimator is the exact C++ local-fraction Hill implementation in
# code/src/local_hill.cpp, loaded here from the project's canonical source.

args <- commandArgs(trailingOnly = TRUE)
output_path <- if (length(args) >= 1L) args[[1L]] else "results/alex_reconciliation.csv"
n_rep <- if (length(args) >= 2L) as.integer(args[[2L]]) else 20L
p <- if (length(args) >= 3L) as.integer(args[[3L]]) else 500L
n <- if (length(args) >= 4L) as.integer(args[[4L]]) else 1000L
stopifnot(n_rep >= 1L, p >= 20L, n >= 100L)

suppressPackageStartupMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)

ar1_gaussian_reconcile <- function(n, p, rho) {
  z <- matrix(NA_real_, n, p)
  z[, 1L] <- rnorm(n)
  innovation_sd <- sqrt(1 - rho^2)
  if (p > 1L) for (j in 2:p)
    z[, j] <- rho * z[, j - 1L] + innovation_sd * rnorm(n)
  z
}

gamma_A <- function(u) 0.5 * exp(-rowSums(u[, 1:4, drop = FALSE]))

generate_y <- function(gamma, uniforms, m = 0.25) {
  y <- exp(gamma * (log(m + 1 - m * uniforms) - log(uniforms)))
  stopifnot(all(is.finite(y)), all(y > 1))
  y
}

# Faithful transcription of Alex lines 420--531.
alex_hill_curve <- function(y, k_values) {
  log_y <- sort(log(y), decreasing = TRUE)
  cumsum(log_y)[k_values] / k_values - log_y[k_values + 1L]
}

alex_select_stable_k <- function(y, k_min = 20L, k_max = floor(length(y) / 2),
                                 window_size = 21L) {
  k_max <- min(as.integer(k_max), length(y) - 2L)
  k_values <- seq.int(as.integer(k_min), k_max)
  hill_values <- alex_hill_curve(y, k_values)
  if (window_size %% 2L == 0L) window_size <- window_size + 1L
  half <- window_size %/% 2L
  stability <- rep(Inf, length(k_values))
  candidates <- seq.int(half + 1L, length(k_values) - half)
  for (ii in candidates) {
    jj <- seq.int(ii - half, ii + half)
    local_k <- k_values[jj]
    local_hill <- hill_values[jj]
    local_mean <- mean(local_hill)
    if (!is.finite(local_mean) || local_mean <= 0) next
    relative_sd <- sd(local_hill) / local_mean
    centered_k <- local_k - mean(local_k)
    slope <- sum(centered_k * local_hill) / sum(centered_k^2)
    relative_slope <- abs(slope) * diff(range(local_k)) / local_mean
    stability[ii] <- relative_sd + relative_slope
  }
  selected <- if (all(!is.finite(stability))) {
    which.min(abs(diff(hill_values))) + 1L
  } else which.min(stability)
  list(k = k_values[selected], hill = hill_values[selected])
}

# Faithful transcription of Alex lines 539--572.
alex_excesses <- function(y, k) {
  ord <- order(y, decreasing = TRUE)
  threshold <- y[ord[k + 1L]]
  exceedance <- rep(FALSE, length(y))
  exceedance[ord[seq_len(k)]] <- TRUE
  log_excess <- numeric(length(y))
  log_excess[exceedance] <- log(y[exceedance] / threshold)
  list(threshold = threshold, exceedance = exceedance, log_excess = log_excess)
}

# Same Epanechnikov calculation as Alex lines 72--248, packaged under a unique
# C++ function name so this file can coexist with other Rcpp sessions.
cppFunction(depends = "Rcpp", includes = '
  #include <algorithm>
  #include <cmath>
  #include <numeric>
  #include <vector>
', code = '
Rcpp::List alex_kernel_scores_cpp(
    const Rcpp::NumericMatrix& z, const Rcpp::NumericVector& log_excess,
    const Rcpp::LogicalVector& exceedance, const double h) {
  const int n = z.nrow(), p = z.ncol();
  Rcpp::NumericVector score(p, NA_REAL), coverage(p, 0.0);
  const int radius = static_cast<int>(std::floor(h * (n + 1.0) + 1e-12));
  const double den = n + 1.0, h2 = h * h;
  std::vector<int> oi(n);
  std::vector<double> r0(n+1), r1(n+1), r2(n+1), w0(n+1), w1(n+1), w2(n+1);
  for (int q=0; q<p; ++q) {
    std::iota(oi.begin(), oi.end(), 0);
    std::sort(oi.begin(), oi.end(), [&z,q](int a,int b){return z(a,q)<z(b,q);});
    r0[0]=r1[0]=r2[0]=w0[0]=w1[0]=w2[0]=0.0;
    for (int j=0; j<n; ++j) {
      int id=oi[j]; double x=(j+1.0)/den, x2=x*x;
      double r=log_excess[id], w=exceedance[id] ? 1.0 : 0.0;
      r0[j+1]=r0[j]+r; r1[j+1]=r1[j]+x*r; r2[j+1]=r2[j]+x2*r;
      w0[j+1]=w0[j]+w; w1[j+1]=w1[j]+x*w; w2[j+1]=w2[j]+x2*w;
    }
    double total=0.0; int valid=0;
    for (int j=0; j<n; ++j) {
      int left=std::max(0,j-radius), upper=std::min(n-1,j+radius)+1;
      double x=(j+1.0)/den, x2=x*x;
      double sr0=r0[upper]-r0[left], sr1=r1[upper]-r1[left], sr2=r2[upper]-r2[left];
      double sw0=w0[upper]-w0[left], sw1=w1[upper]-w1[left], sw2=w2[upper]-w2[left];
      double num=sr0-(sr2-2*x*sr1+x2*sr0)/h2;
      double dst=sw0-(sw2-2*x*sw1+x2*sw0)/h2;
      if (R_finite(num) && R_finite(dst) && dst>1e-12) {total += num/dst; ++valid;}
    }
    if (valid>0) {score[q]=total/valid; coverage[q]=valid/(double)n;}
  }
  return Rcpp::List::create(Rcpp::Named("score")=score,
                            Rcpp::Named("coverage")=coverage);
}')

article_scores <- function(z, y, h, alpha, epsilon = 0.05) {
  vapply(seq_len(ncol(z)), function(j)
    unname(score_coordinate_cpp(z[, j], y, h, alpha, epsilon)[["score"]]),
    numeric(1L))
}

alex_scores <- function(z, y, h, selected_k = NULL) {
  if (is.null(selected_k)) selected_k <- alex_select_stable_k(y)$k
  excess <- alex_excesses(y, selected_k)
  fit <- alex_kernel_scores_cpp(z, excess$log_excess, excess$exceedance, h)
  list(scores = fit$score, coverage = fit$coverage, k = selected_k,
       threshold = excess$threshold)
}

# Independent direct formula used only as a numerical audit of the optimized
# Alex C++ cumulative-sum formula. It follows the mathematical kernel ratio,
# not the implementation's prefix-sum algebra.
alex_score_direct <- function(x, log_excess, exceedance, h) {
  ord <- order(x)
  r <- log_excess[ord]
  w <- as.numeric(exceedance[ord])
  grid <- seq_along(x) / (length(x) + 1)
  values <- vapply(grid, function(u) {
    weights <- pmax(0, 1 - ((grid - u) / h)^2)
    denominator <- sum(weights * w)
    if (denominator <= 1e-12) return(NA_real_)
    sum(weights * r) / denominator
  }, numeric(1L))
  mean(values[is.finite(values)])
}

selection_metrics <- function(scores, d) {
  ordering <- order(scores, seq_along(scores), na.last = TRUE)
  ranks <- match(1:4, ordering)
  selected <- ordering[seq_len(d)]
  c(sure = as.integer(all(1:4 %in% selected)),
    exact = as.integer(d == 4L && setequal(selected, 1:4)),
    max_active_rank = max(ranks), mean_active_rank = mean(ranks))
}

rho <- 0.25
m <- 0.25
alpha <- n^(-0.30)
h_article <- n^(-0.20) / 2
h_alex <- 0.20
base_seed <- 731921L
rows <- vector("list", n_rep * 2L * 4L * 2L)
at <- 0L
audit_max_abs <- numeric(n_rep * 2L)
audit_at <- 0L

started <- proc.time()[[3L]]
for (replicate in seq_len(n_rep)) {
  # Each truth-DGP variant uses precisely the same latent AR(1) matrix and the
  # same response uniforms within a replication.
  set.seed(base_seed + 1009L * replicate)
  z <- ar1_gaussian_reconcile(n, p, rho)
  response_uniforms <- runif(n)
  u_population <- pnorm(z[, 1:4, drop = FALSE])
  u_ranks <- apply(z[, 1:4, drop = FALSE], 2L, rank, ties.method = "first") / (n + 1)

  for (truth_dgp in c("pnorm", "empirical_rank")) {
    u_truth <- if (truth_dgp == "pnorm") u_population else u_ranks
    y <- generate_y(gamma_A(u_truth), response_uniforms, m)
    selected_k <- alex_select_stable_k(y)$k

    fits <- list(
      article_article_h = list(scores = article_scores(z, y, h_article, alpha),
                               h = h_article, k = NA_integer_),
      article_alex_h = list(scores = article_scores(z, y, h_alex, alpha),
                            h = h_alex, k = NA_integer_),
      alex_global_article_h = c(alex_scores(z, y, h_article, selected_k),
                                list(h = h_article)),
      alex_global_kernel = c(alex_scores(z, y, h_alex, selected_k), list(h = h_alex))
    )

    # Verify the optimized Alex formula on the first coordinate for every DGP
    # and replication. This comparison also detects boundary/radius mistakes.
    ex <- alex_excesses(y, selected_k)
    direct <- alex_score_direct(z[, 1L], ex$log_excess, ex$exceedance, h_alex)
    audit_at <- audit_at + 1L
    audit_max_abs[audit_at] <- abs(direct - fits$alex_global_kernel$scores[1L])

    for (estimator in names(fits)) for (d in c(4L, 20L)) {
      fit <- fits[[estimator]]
      metrics <- selection_metrics(fit$scores, d)
      at <- at + 1L
      rows[[at]] <- data.frame(
        replicate = replicate, seed = base_seed + 1009L * replicate,
        n = n, p = p, rho = rho, truth_dgp = truth_dgp,
        estimator = estimator, h = fit$h, alpha = alpha,
        selected_k = fit$k, d = d,
        sure = metrics[["sure"]], exact = metrics[["exact"]],
        max_active_rank = metrics[["max_active_rank"]],
        mean_active_rank = metrics[["mean_active_rank"]],
        min_score = min(fit$scores, na.rm = TRUE),
        max_score = max(fit$scores, na.rm = TRUE),
        stringsAsFactors = FALSE)
    }
  }
  if (replicate == 1L || replicate %% 5L == 0L || replicate == n_rep)
    message(sprintf("reconciliation replication %d/%d", replicate, n_rep))
}

out <- do.call(rbind, rows[seq_len(at)])
out$formula_audit_max_abs <- max(audit_max_abs[seq_len(audit_at)])
out$elapsed_total_seconds <- proc.time()[[3L]] - started
stopifnot(out$formula_audit_max_abs[[1L]] < 1e-10)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(out, output_path, row.names = FALSE)
message(sprintf("WROTE %s (%d rows); formula max abs error %.3g",
                output_path, nrow(out), out$formula_audit_max_abs[[1L]]))
