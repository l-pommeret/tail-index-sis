## The fast grid scorer must reproduce the reference Yoshida--Umezu
## implementation exactly, for every (h, k) cell.  Only then may the grid
## search be trusted to select a tuning on the competitor's behalf.
source("code/R/yoshida_umezu.R")
source("code/R/yoshida_umezu_grid.R")
source("code/R/generate5.R")

fail <- 0L
report <- function(ok, what) {
  cat(if (ok) "PASS " else "FAIL ", what, "\n", sep = "")
  if (!ok) fail <<- fail + 1L
}

set.seed(20260817)
EVAL <- seq(.02, .98, length.out = 25)

for (spec in list(list(n = 600L, p = 40L, model = "A1"),
                  list(n = 800L, p = 25L, model = "A2"),
                  list(n = 600L, p = 30L, model = "B1"))) {
  d <- simulate_dataset5(spec$n, spec$p, 0.25, spec$model, 4242L + spec$n)
  u <- apply(d$z, 2L, function(x)
    rank(x, ties.method = "average") / (spec$n + 1))
  ks <- as.integer(round(c(.05, .06, .072, .08) * spec$n))
  for (h in c(0.5, 1, 1.5, 2, 3)) {
    fast <- yu_score_matrix_grid(u, d$y, ks, h, EVAL)
    for (ik in seq_along(ks)) {
      ref <- yu_score_matrix(u, d$y, ks[ik], h, EVAL)
      same_s <- identical(all.equal(ref$scores, fast$scores[, ik],
                                    tolerance = 0), TRUE)
      same_u <- identical(all.equal(ref$undefined, fast$undefined[, ik],
                                    tolerance = 0), TRUE)
      report(same_s && same_u,
             sprintf("%s n=%d h=%.1f k=%d : grid == reference",
                     spec$model, spec$n, h, ks[ik]))
    }
  }
}

## The selection rule must be a pure function of the tuning table.
tab <- data.frame(k = c(100L, 100L, 120L, 120L), h = c(1, 2, 1, 2),
                  sure20 = c(0.4, 0.7, 0.7, 0.2), ermax = c(50, 30, 20, 90))
source("code/R/yu_select.R")
best <- yu_select_best(tab)
report(best$k == 120L && best$h == 1,
       "selection takes the Sure-20 argmax, ties broken by E(Rmax)")

tab2 <- data.frame(k = c(100L, 120L), h = c(1, 2),
                   sure20 = c(0.7, 0.7), ermax = c(40, 25))
best2 <- yu_select_best(tab2)
report(best2$k == 120L && best2$h == 2, "tie on Sure-20 broken by E(Rmax)")

if (fail) stop(fail, " YU-grid check(s) failed") else
  cat("ALL YU-GRID CHECKS PASSED\n")
