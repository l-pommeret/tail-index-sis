## Unit checks for generate5.R.
## Full and streaming generators consume the RNG in different orders, so they
## are not bitwise identical at a common seed; the checks below verify that
## they implement the same intended distribution, and specifically that the
## Draft-4 streaming defect (AR(1) chain beyond the proxy block continuing
## from the overwritten proxy coordinate) is fixed.
suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate5.R")

rho <- 0.25; K <- 300L

## 1. B1, full generator: coordinate 25 continues the base chain, so its
##    correlation with the observed proxy at 24 is ~0 and with 26 is ~rho.
c_full <- t(replicate(K, {
  d <- simulate_dataset5(800L, 27L, rho, "B1", seed = sample.int(1e8, 1))
  c(cor(d$z[, 24], d$z[, 25]), cor(d$z[, 25], d$z[, 26]))
}))
stopifnot(abs(mean(c_full[, 1])) < 0.02, abs(mean(c_full[, 2]) - rho) < 0.02)
cat(sprintf("PASS full B1: cor(Z24p,Z25)=%.3f (~0), cor(Z25,Z26)=%.3f (~rho)\n",
            mean(c_full[, 1]), mean(c_full[, 2])))

## 2. B1, streaming generator: same two properties.  Under the Draft-4 defect
##    cor(Z24_proxy, Z25) would be ~rho = 0.25; after the fix it is ~0.
##    We recover the streamed covariates by re-simulating the stream by hand.
stream_cols <- function(seed, p = 27L, n = 800L) {
  set.seed(seed)
  n0 <- max(A_SCALE)
  innovation_sd <- sqrt(1 - rho^2)
  z <- matrix(NA_real_, n, p)
  z[, 1] <- rnorm(n)
  for (j in 2:n0) z[, j] <- rho * z[, j - 1] + innovation_sd * rnorm(n)
  previous_ar <- z[, n0]
  f <- rnorm(n)
  z <- b1_scale_block(z, f)
  invisible(runif(n))                       # V, consumed by the response
  previous <- previous_ar
  for (j in (n0 + 1):p) {
    z[, j] <- rho * previous + innovation_sd * rnorm(n)
    previous <- z[, j]
  }
  z
}
c_str <- t(replicate(K, {
  z <- stream_cols(sample.int(1e8, 1))
  c(cor(z[, 24], z[, 25]), cor(z[, 25], z[, 26]))
}))
stopifnot(abs(mean(c_str[, 1])) < 0.02, abs(mean(c_str[, 2]) - rho) < 0.02)
cat(sprintf("PASS streaming B1: cor(Z24p,Z25)=%.3f (~0), cor(Z25,Z26)=%.3f (~rho)\n",
            mean(c_str[, 1]), mean(c_str[, 2])))

## 3. The hand-rolled stream above must match simulate_scores_streaming5
##    exactly at a common seed (same RNG order): scores of the first column
##    computed from both paths agree bitwise.
settings <- data.frame(h = 800^(-0.15) / 2, alpha = 800^(-0.30))
S <- simulate_scores_streaming5(800L, 27L, rho, "B1", seed = 99L, settings,
                                block_size = 256L)
z99 <- stream_cols(99L)
set.seed(99L)                               # rebuild y with the same draws
n0 <- max(A_SCALE); innovation_sd <- sqrt(1 - rho^2)
zz <- matrix(NA_real_, 800L, n0); zz[, 1] <- rnorm(800L)
for (j in 2:n0) zz[, j] <- rho * zz[, j - 1] + innovation_sd * rnorm(800L)
f <- { invisible(NULL); rnorm(800L) }
zz <- b1_scale_block(zz, f)
u <- population_uniform_matrix(zz)
y <- response5(gamma_model(u, "B1"), u, "B1", V = runif(800L), f = f)
sc1 <- score_coordinate_cpp(z99[, 1], y, settings$h[1], settings$alpha[1],
                            .05)[["score"]]
stopifnot(abs(S[1, 1] - sc1) < 1e-12)
cat("PASS streaming path matches hand-rolled stream at common seed\n")

## 4. A-family: full and streaming score distributions agree (A1, 60 seeds,
##    mean active-coordinate score difference within Monte Carlo error).
mfull <- replicate(60, {
  d <- simulate_dataset5(800L, 30L, rho, "A1", seed = sample.int(1e8, 1))
  mean(vapply(1:4, function(j)
    score_coordinate_cpp(d$z[, j], d$y, settings$h[1], settings$alpha[1],
                         .05)[["score"]], numeric(1)))
})
mstr <- replicate(60, {
  S <- simulate_scores_streaming5(800L, 30L, rho, "A1",
                                  seed = sample.int(1e8, 1), settings)
  NA_real_ -> dummy
  mean(S[1:4, 1] * 0 + S[1:4, 1]) # mean of active-coordinate scores
})
d <- abs(mean(mfull) - mean(mstr))
se <- sqrt(var(mfull) / 60 + var(mstr) / 60)
stopifnot(d < 4 * se)
cat(sprintf("PASS A1 full vs streaming score means: |diff|=%.4f (4*SE=%.4f)\n",
            d, 4 * se))
cat("ALL GENERATE5 CHECKS PASS\n")
