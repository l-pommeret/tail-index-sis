## Grid search at rho = 0 over (n, p), all four models, comparing:
##   1 screen propose au reglage publie
##   2 screen agrege par rang median sur la grille 3x3
##   3 screen agrege par rang minimum sur la grille 3x3
##   4 quantile SIS tau = 0.95
##   5 pipeline : quantile SIS top-25 puis rang minimum sur la grille 3x3
##
## Fresh data (n = 4000 is outside the campaign grid), seed stream
## 81000019 + cell*100003 + arm*10007 + r*307, disjoint from every other stream.
## Cells are run longest-first, each as one mclapply over model x replicate so
## every worker stays busy, and each cell is checkpointed on completion.
##
## usage: KAPPA=0.20 Rscript code/R/grid_search_r0.R OUTDIR [NREP] [CORES]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R"); source("code/R/qa_sis.R")

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[1L] else "results/grid_r0"
NREP  <- if (length(args) >= 2L) as.integer(args[2L]) else 40L
CORES <- if (length(args) >= 3L) as.integer(args[3L]) else 92L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.05, 0.10, 0.15)
RHO <- 0
D1 <- 25L
MODELS <- c("M1", "M2", "M3", "M4")
NS <- c(1000L, 2000L, 4000L)
PS <- c(200L, 500L, 1000L, 2000L)
RULES <- c("screen seul", "screen 9 mediane", "screen 9 minimum",
           "quantile .95", "pipeline + 9 min")

## per-replicate cost proxy (seconds), for longest-first scheduling
cost1 <- c("1000" = 3.6, "2000" = 7.5, "4000" = 29)

cells <- expand.grid(p = PS, n = NS, KEEP.OUT.ATTRS = FALSE)
cells$id <- sprintf("r0_n%04d_p%04d", cells$n, cells$p)
cells$cost <- cost1[as.character(cells$n)] * cells$p / 1000 * NREP * length(MODELS)
cells <- cells[order(-cells$cost), ]

score_vec <- function(z, y, n, p, a, b) {
  h <- n^(-b) / 2; alpha <- n^(-a)
  vapply(seq_len(p), function(j)
    score_coordinate_cpp(z[, j], y, h, alpha, EPS), numeric(5))["score", ]
}

run_one <- function(n, p, model, r, cellseed) {
  seed <- cellseed + match(model, MODELS) * 10007L + r * 307L
  d <- simulate_dataset3(n, p, RHO, model, seed)
  base <- score_vec(d$z, d$y, n, p, ASTAR, BSTAR)
  R <- matrix(NA_real_, p, length(AGRID) * length(BGRID)); ii <- 0L
  for (a in AGRID) for (b in BGRID) {
    ii <- ii + 1L
    R[, ii] <- rank(score_vec(d$z, d$y, n, p, a, b), ties.method = "first",
                    na.last = TRUE)
  }
  amed <- apply(R, 1L, median); amin <- apply(R, 1L, min)
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (n + 1))
  ord_q <- order(-qa_sis_scores(uh, d$y, tau = 0.95), seq_len(p))
  surv <- ord_q[seq_len(min(D1, p))]
  r_pipe <- if (!all(1:4 %in% surv)) {
    Inf
  } else {
    max(match(1:4, surv[order(amin[surv], surv)]))
  }
  list(model = model, n = n, p = p, rho = RHO, replicate = r, seed = seed,
       rmax = c(max(match(1:4, order(base, seq_len(p), na.last = TRUE))),
                max(match(1:4, order(amed, seq_len(p)))),
                max(match(1:4, order(amin, seq_len(p)))),
                max(match(1:4, ord_q)),
                r_pipe))
}

cat(sprintf("%d cellules, %.1f core-heures estimees\n", nrow(cells),
            sum(cells$cost) / 3600)); flush.console()
for (i in seq_len(nrow(cells))) {
  path <- file.path(outdir, paste0(cells$id[i], ".rds"))
  if (file.exists(path)) next
  n <- cells$n[i]; p <- cells$p[i]
  cellseed <- 81000019L + i * 100003L
  grid <- expand.grid(model = MODELS, r = seq_len(NREP),
                      stringsAsFactors = FALSE)
  out <- parallel::mclapply(seq_len(nrow(grid)), function(k)
    run_one(n, p, grid$model[k], grid$r[k], cellseed),
    mc.cores = CORES, mc.preschedule = FALSE)
  bad <- vapply(out, function(z) inherits(z, "try-error") || is.null(z),
                logical(1))
  if (any(bad)) stop(cells$id[i], ": ", sum(bad), " failed: ",
                     as.character(out[bad][[1]]))
  saveRDS(out, path, compress = "xz")
  cat(format(Sys.time(), "%H:%M:%S"), "cellule", cells$id[i], "faite\n")
  flush.console()
}

## ------------------------------------------------------------------ tables --
fm <- function(v) formatC(v, format = "f", digits = 3)
files <- sort(list.files(outdir, pattern = "^r0_.*\\.rds$", full.names = TRUE))
S <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f)
  do.call(rbind, lapply(MODELS, function(m) {
    w <- x[vapply(x, function(z) z$model == m, logical(1))]
    do.call(rbind, lapply(seq_along(RULES), function(k) {
      v <- vapply(w, function(z) z$rmax[k], numeric(1))
      data.frame(n = w[[1]]$n, p = w[[1]]$p, model = m, rule = RULES[k],
                 sure4 = mean(v <= 4), sure20 = mean(v <= 20))
    }))
  }))
}))
write.csv(S, file.path(outdir, "summary.csv"), row.names = FALSE)

cat("\n=== Sure-4 / Sure-20, moyenne sur les 4 modeles, rho = 0 ===\n\n")
for (nn in NS) {
  cat(sprintf("n = %d\n", nn))
  cat(sprintf("%-20s %s\n", "regle",
              paste(sprintf("%-16s", paste0("p=", PS)), collapse = "")))
  for (k in RULES) {
    v <- vapply(PS, function(pp) {
      w <- S[S$n == nn & S$p == pp & S$rule == k, ]
      if (!nrow(w)) return("  -")
      sprintf("%s/%s", fm(mean(w$sure4)), fm(mean(w$sure20)))
    }, character(1))
    cat(sprintf("%-20s %s\n", k, paste(sprintf("%-16s", v), collapse = "")))
  }
  cat("\n")
}
