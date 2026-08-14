options(warn = 2)
source("code/R/estimator.R")
source("code/R/generate.R")
source("code/R/streaming.R")
library(Rcpp)
sourceCpp("code/src/local_hill.cpp")

set.seed(9127)
n <- 180L
u <- rank(runif(n)) / (n + 1)
y <- exp(rexp(n, rate = 2))
ev <- u[u >= .05 & u <= .95]
r <- local_hill_r(u, y, ev, h = .1, alpha = .2)
c <- local_hill_cpp(u, y, ev, h = .1, alpha = .2)
stopifnot(isTRUE(all.equal(r$estimate, c$estimate, tolerance = 1e-12)),
          identical(r$local_n, c$local_n), identical(r$local_k, c$local_k),
          identical(r$underpopulated, c$underpopulated))
sr <- score_coordinate_r(u, y, h = .1, alpha = .2, epsilon = .05)
sc <- score_coordinate_cpp(u, y, h = .1, alpha = .2, epsilon = .05)
stopifnot(isTRUE(all.equal(unname(sr), unname(sc), tolerance = 1e-12)))

# Underpopulation is never encoded as the best possible (zero) score.
tiny_r <- score_coordinate_r(1:8, exp(1:8), h = .01, alpha = .1, epsilon = .05)
tiny_c <- score_coordinate_cpp(1:8, exp(1:8), h = .01, alpha = .1, epsilon = .05)
stopifnot(is.infinite(tiny_r[["score"]]), is.infinite(tiny_c[["score"]]),
          tiny_r[["n_eval"]] == 0, tiny_c[["n_eval"]] == 0)

d0 <- simulate_dataset(100, 20, 0, "N", seed = 17, block_size = 3)
stopifnot(all(d0$y > 1), all(d0$gamma == .3),
          all(d0$u > 0), all(d0$u < 1),
          identical(d0$u, pnorm(d0$x)))
d1 <- simulate_dataset(100, 20, .75, "A", seed = 18, block_size = 7)
stopifnot(all(d1$y > 1), all(d1$gamma > 0), ncol(d1$x) == 20)
dF <- simulate_dataset(100, 20, .5, "F", seed = 19, block_size = 11)
stopifnot(all(dF$y > 1), isTRUE(all.equal(dF$gamma,
  gamma_model(dF$u, "A"), tolerance = 0)))

# AR(1) recursion is invariant to the declared generation block size.
set.seed(812); za <- ar1_gaussian(50, 31, .5, block_size = 1)
set.seed(812); zb <- ar1_gaussian(50, 31, .5, block_size = 13)
stopifnot(identical(za, zb))

# Full simulation and scores are invariant to processing block size.
ss1 <- simulate_score_streaming(120, 23, .5, "F", 1, 7321, .1, .2, .05, 1)
ss2 <- simulate_score_streaming(120, 23, .5, "F", 1, 7321, .1, .2, .05, 11)
stopifnot(identical(ss1, ss2), all(ss1$scores > 0), ss1$y_min > 1)

# A population PIT is not the deterministic empirical-rank grid.
set.seed(115)
zp <- matrix(rnorm(800), 200, 4)
up <- population_uniform_matrix(zp)
ur <- rank_uniform_matrix(zp)
stopifnot(identical(up, pnorm(zp)), !isTRUE(all.equal(up, ur, tolerance = 0)))

# Model-level mathematical invariants.
grid <- as.matrix(expand.grid(rep(list(seq(0, 1, length.out = 7)), 4)))
stopifnot(all(gamma_model(grid, "C") > 0), all(gamma_model(grid, "D") > 0),
          all(gamma_model(grid, "E") > 0))
cat("PASS core estimator and generator tests\n")
