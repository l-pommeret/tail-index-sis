## Timing probe for the extended comparison grid.
## usage: Rscript code/R/bench_grid.R
suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R")
source("code/R/qa_sis.R")
source("code/R/yoshida_umezu.R")
EPS <- 0.05; ASTAR <- 0.30; BSTAR <- 0.10
TAUS <- c(0.90, 0.95, 0.975, 0.99)

probe <- function(n, p, rho = 0.25, m = "M1") {
  tt <- function(expr) unname(system.time(expr)[3])
  t_gen <- tt(d <- simulate_dataset3(n, p, rho, m, 12345L))
  t_pit <- tt(uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average")/(n+1)))
  h <- n^(-BSTAR)/2; alpha <- n^(-ASTAR)
  t_ours <- tt(sc <- vapply(seq_len(p), function(j)
    score_coordinate_cpp(d$z[,j], d$y, h, alpha, EPS), numeric(5))["score",])
  t_yu <- tt(yu <- yu_score_matrix(uh, d$y, k = floor(0.072*n), h = 1))
  t_q <- tt(q1 <- qa_sis_scores(uh, d$y, tau = 0.95))
  cat(sprintf("n=%5d p=%5d | gen %6.2f pit %6.2f ours %6.2f yu %7.2f q(1tau) %7.2f -> total(4taus) %7.2f s\n",
              n, p, t_gen, t_pit, t_ours, t_yu, t_q,
              t_gen + t_pit + t_ours + t_yu + 4*t_q))
  invisible(NULL)
}

for (n in c(1000L, 2000L, 5000L)) for (p in c(200L, 1000L)) probe(n, p)
probe(5000L, 2000L)
