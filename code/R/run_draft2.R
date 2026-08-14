## Draft-2 experiment driver.
## usage: Rscript code/R/run_draft2.R {tuning|sensitivity|comparison} OUTPUT.rds CORES
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
what <- args[1L]; output_path <- args[2L]; cores <- as.integer(args[3L])

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate2.R")

EPS <- 0.05; BLOCK <- 256L

run_jobs <- function(jobs, fn) {
  out <- parallel::mclapply(seq_len(nrow(jobs)), function(i) fn(jobs[i, ]),
                            mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(out, inherits, logical(1), "try-error")
  if (any(bad)) stop("failed jobs: ", sum(bad))
  out
}

if (what == "tuning") {
  grid <- expand.grid(model = c("M1", "M2", "M3", "M4"),
                      a = c(.25, .30, .35, .40, .45, .50),
                      b = c(0, .05, .10, .15, .20, .30, .40),
                      stringsAsFactors = FALSE)
  jobs <- grid[rep(seq_len(nrow(grid)), each = 200L), ]
  jobs$replicate <- rep(1:200, nrow(grid))
  jobs$seed <- 2000003L + as.integer(interaction(grid, drop = TRUE))[
    rep(seq_len(nrow(grid)), each = 200L)] * 100003L + jobs$replicate * 211L
  fn <- function(z) {
    n <- 2000L; p <- 1000L; rho <- 0.25
    h <- n^(-z$b) / 2; alpha <- n^(-z$a)
    d <- simulate_score_streaming2(n, p, rho, z$model, z$seed, h, alpha, EPS, BLOCK)
    o <- order(d$scores, seq_along(d$scores), na.last = TRUE)
    list(model = z$model, a = z$a, b = z$b, replicate = z$replicate,
         seed = z$seed, active_ranks = match(1:4, o))
  }
} else if (what == "sensitivity") {
  cells <- rbind(
    data.frame(n = c(1000L, 2000L, 5000L), p = 1000L, rho = 0.25),
    data.frame(n = 2000L, p = c(100L, 500L, 5000L), rho = 0.25),
    data.frame(n = 2000L, p = 1000L, rho = c(0, 0.5, 0.75)))
  jobs <- cells[rep(seq_len(nrow(cells)), each = 200L), ]
  jobs$replicate <- rep(1:200, nrow(cells))
  jobs$seed <- 3000017L + rep(seq_len(nrow(cells)), each = 200L) * 100019L +
    jobs$replicate * 307L
  ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
  BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
  fn <- function(z) {
    h <- z$n^(-BSTAR) / 2; alpha <- z$n^(-ASTAR)
    d <- simulate_score_streaming2(z$n, z$p, z$rho, "M1", z$seed, h, alpha, EPS, BLOCK)
    o <- order(d$scores, seq_along(d$scores), na.last = TRUE)
    list(n = z$n, p = z$p, rho = z$rho, replicate = z$replicate, seed = z$seed,
         a = ASTAR, b = BSTAR, active_ranks = match(1:4, o), ordering = o[1:500])
  }
} else if (what == "comparison") {
  source("code/R/qa_sis.R"); source("code/R/yoshida_umezu.R")
  grid <- data.frame(model = c("M1", "M2", "M3", "M4"))
  jobs <- grid[rep(1:4, each = 200L), , drop = FALSE]
  jobs$replicate <- rep(1:200, 4L)
  jobs$seed <- 4000037L + rep(1:4, each = 200L) * 100043L + jobs$replicate * 401L
  ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
  BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
  fn <- function(z) {
    n <- 2000L; p <- 1000L; rho <- 0.25
    d <- simulate_dataset2(n, p, rho, z$model, z$seed, BLOCK)
    uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (n + 1))
    ## proposed method (increasing scores)
    h <- n^(-BSTAR) / 2; alpha <- n^(-ASTAR)
    sc <- vapply(seq_len(p), function(j)
      score_coordinate_cpp(d$z[, j], d$y, h, alpha, EPS), numeric(5))["score", ]
    o_ours <- order(sc, seq_len(p), na.last = TRUE)
    ## Yoshida-Umezu (decreasing utilities), their baseline tuning
    yu <- yu_score_matrix(uh, d$y, k = floor(0.072 * n), h = 1)
    o_yu <- order(-yu$scores, seq_len(p))
    ## quantile-adaptive SIS at tau = .95 (decreasing utilities)
    qs95 <- qa_sis_scores(uh, d$y, tau = 0.95)
    o_q95 <- order(-qs95, seq_len(p))
    ## supplementary: tau = .99
    qs99 <- qa_sis_scores(uh, d$y, tau = 0.99)
    o_q99 <- order(-qs99, seq_len(p))
    rk <- function(o) match(1:4, o)
    list(model = z$model, replicate = z$replicate, seed = z$seed,
         ranks_ours = rk(o_ours), ranks_yu = rk(o_yu),
         ranks_q95 = rk(o_q95), ranks_q99 = rk(o_q99),
         top20_q95 = o_q95[1:20], top20_q99 = o_q99[1:20])
  }
} else stop("unknown experiment")

meta <- list(what = what, started_utc = format(Sys.time(), tz = "UTC"),
             astar = Sys.getenv("ASTAR", ""), bstar = Sys.getenv("BSTAR", ""),
             r_version = R.version.string, cores = cores)
out <- run_jobs(jobs, fn)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(meta = meta, jobs = out), output_path, compress = "xz")
cat("WROTE", output_path, "JOBS", length(out), "\n")
