## Where does the comparison time go, and is the quantile-SIS LP solver
## replaceable?  Barrodale--Roberts (rq.fit.br, used by code/R/qa_sis.R) versus
## Frisch--Newton (rq.fit.fnb), on the same B-spline design.
suppressPackageStartupMessages({ library(quantreg); library(splines) })
set.seed(4321)
for (n in c(1000L, 2000L, 5000L)) {
  x <- rnorm(n); y <- exp(rnorm(n))
  B <- cbind(1, bs(x, df = 3))
  tt <- function(e) unname(system.time(e)[3])
  m <- 40L
  t_br  <- tt(for (i in 1:m) fb <- rq.fit.br(B, y, tau = 0.95))
  t_fn  <- tt(for (i in 1:m) ff <- rq.fit.fnb(B, y, tau = 0.95))
  q <- as.numeric(quantile(y, 0.95, type = 7))
  w_br <- mean((B %*% fb$coefficients - q)^2)
  w_fn <- mean((B %*% ff$coefficients - q)^2)
  cat(sprintf("n=%5d  br %6.1f ms/fit   fnb %6.1f ms/fit   speedup %4.1fx",
              n, 1000 * t_br / m, 1000 * t_fn / m, t_br / t_fn),
      sprintf("   utility rel.diff %.2e\n", abs(w_br - w_fn) / w_br))
}
