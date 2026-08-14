## Campagne finale sur la suite de modèles figée (code/R/generate4.R).
##
##   n   = 2000
##   p   in {500, 1000, 2000}
##   rho = 0.25
##   modèles M1, M2, M3, M4
##   1000 réplications par cellule et par modèle
##
## Règles comparées, toutes sur les mêmes jeux de données :
##   1 screen proposé au réglage retenu (a*, b*) = (0.30, 0.15), choisi par la
##     phase de réglage code/R/tuning4.R sur cette même suite de modèles
##   2 screen agrégé par rang minimum sur la grille 3x3
##   3 Yoshida--Umezu, réglage publié h = 1, k = floor(0.072 n)
##   4-6 quantile SIS a tau = 0.90, 0.95, 0.99
##
## Le réglage retenu appartenant à la grille 3x3, la règle 1 réutilise cette
## passe du score au lieu d'en recalculer une.
##
## Métriques : Sure-d pour d dans {4,10,20,30,50}, E(Rmax) et Med(Rmax).
##
## Flux de graines 131000021 + cellule*100003 + modele*10007 + r*307, disjoint
## de tous les autres.
##
## usage: Rscript code/R/campaign5.R OUTDIR [NREP] [CORES]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate4.R"); source("code/R/qa_sis.R")
source("code/R/yoshida_umezu.R")

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[1L] else "results/campaign5"
NREP  <- if (length(args) >= 2L) as.integer(args[2L]) else 1000L
CORES <- if (length(args) >= 3L) as.integer(args[3L]) else 92L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.15"))
## grille d'agregation : bloc 3x3 de bons reglages adjacent a l'optimum.
## a = 0.25 en est exclu deliberement, la ligne s'effondrant pour M1, M3 et M4
## (Sure-20 au plus 0.110, 0.085 et 0.250) ; y inclure des reglages mauvais
## degrade l'agregation au lieu de l'ameliorer.
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.10, 0.15, 0.20)
N <- 2000L
PS <- c(500L, 1000L, 2000L)
RHO <- 0.25
MODELS <- c("M1", "M2", "M3", "M4")
TAUS <- c(0.90, 0.95, 0.99)
DS <- c(4, 10, 20, 30, 50)
RULES <- c("screen seul", "screen 9 minimum", "Yoshida-Umezu",
           paste0("quantile .", c("90", "95", "99")))
## indice du réglage publié dans la grille, parcourue a-externe / b-interne
IDX_PUB <- which(rep(AGRID, each = length(BGRID)) == ASTAR &
                 rep(BGRID, times = length(AGRID)) == BSTAR)
stopifnot(length(IDX_PUB) == 1L)

cells <- data.frame(p = PS)
cells$id <- sprintf("c5_n%04d_p%04d_r%03.0f", N, cells$p, RHO * 100)
cells$cost <- 16.6 * cells$p / 1000 * NREP * length(MODELS)
cells <- cells[order(-cells$cost), ]

score_vec <- function(z, y, p, a, b) {
  h <- N^(-b) / 2; alpha <- N^(-a)
  vapply(seq_len(p), function(j)
    score_coordinate_cpp(z[, j], y, h, alpha, EPS), numeric(5))["score", ]
}

run_one <- function(p, model, r, cellseed) {
  seed <- cellseed + match(model, MODELS) * 10007L + r * 307L
  d <- simulate_dataset4(N, p, RHO, model, seed)
  S <- matrix(NA_real_, p, length(AGRID) * length(BGRID))
  R <- S; ii <- 0L
  for (a in AGRID) for (b in BGRID) {
    ii <- ii + 1L
    S[, ii] <- score_vec(d$z, d$y, p, a, b)
    R[, ii] <- rank(S[, ii], ties.method = "first", na.last = TRUE)
  }
  base <- S[, IDX_PUB]                       # réglage publié, déjà calculé
  amin <- apply(R, 1L, min)
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (N + 1))
  yu <- yu_score_matrix(uh, d$y, k = floor(0.072 * N), h = 1)
  ords <- c(list(order(base, seq_len(p), na.last = TRUE),
                 order(amin, seq_len(p)),
                 order(-yu$scores, seq_len(p))),
            lapply(TAUS, function(tt)
              order(-qa_sis_scores(uh, d$y, tau = tt), seq_len(p))))
  list(model = model, n = N, p = p, rho = RHO, replicate = r, seed = seed,
       rmax = vapply(ords, function(o) max(match(A_GAMMA, o)), numeric(1)),
       yu_undefined = mean(!is.finite(yu$scores)),
       t4g = vapply(ords, function(o) sum(o[1:4] %in% A_GAMMA), numeric(1)),
       t4s = vapply(ords, function(o) sum(o[1:4] %in% A_SCALE), numeric(1)),
       t24s = vapply(ords, function(o) sum(o[1:24] %in% A_SCALE), numeric(1)))
}

cat(sprintf("%d cellules x %d modeles x %d replications, %.1f core-heures estimees\n",
            nrow(cells), length(MODELS), NREP, sum(cells$cost) / 3600))
flush.console()
for (i in seq_len(nrow(cells))) {
  path <- file.path(outdir, paste0(cells$id[i], ".rds"))
  if (file.exists(path)) next
  cellseed <- 131000021L + i * 100003L
  grid <- expand.grid(model = MODELS, r = seq_len(NREP),
                      stringsAsFactors = FALSE)
  out <- parallel::mclapply(seq_len(nrow(grid)), function(k)
    run_one(cells$p[i], grid$model[k], grid$r[k], cellseed),
    mc.cores = CORES, mc.preschedule = FALSE)
  bad <- vapply(out, function(z) inherits(z, "try-error") || is.null(z),
                logical(1))
  if (any(bad)) stop(cells$id[i], ": ", sum(bad), " echecs : ",
                     as.character(out[bad][[1]]))
  saveRDS(out, path, compress = "xz")
  cat(format(Sys.time(), "%H:%M:%S"), "cellule", cells$id[i], "faite\n")
  flush.console()
}

## ------------------------------------------------------------------ tables --
files <- sort(list.files(outdir, pattern = "^c5_.*\\.rds$", full.names = TRUE))
S <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f)
  do.call(rbind, lapply(MODELS, function(m) {
    w <- x[vapply(x, function(z) z$model == m, logical(1))]
    do.call(rbind, lapply(seq_along(RULES), function(k) {
      v <- vapply(w, function(z) z$rmax[k], numeric(1))
      r <- data.frame(n = N, p = w[[1]]$p, rho = RHO, model = m,
                      rule = RULES[k], reps = length(w))
      for (d in DS) r[[paste0("sure", d)]] <- mean(v <= d)
      r$ermax <- mean(v)
      r$medmax <- median(v)
      r$top4_gamma <- mean(vapply(w, function(z) z$t4g[k], numeric(1)))
      r$top4_scale <- mean(vapply(w, function(z) z$t4s[k], numeric(1)))
      r$top24_scale <- mean(vapply(w, function(z) z$t24s[k], numeric(1)))
      r
    }))
  }))
}))
write.csv(S, file.path(outdir, "summary.csv"), row.names = FALSE)

fm <- function(v, k = 3) formatC(v, format = "f", digits = k)
cat(sprintf("\n=== n=%d, rho=%.2f, %d replications (erreur type <= %.4f) ===\n",
            N, RHO, NREP, 0.5 / sqrt(NREP)))
for (m in MODELS) {
  cat(sprintf("\n%s\n", m))
  cat(sprintf("%-20s %s\n", "regle",
              paste(sprintf("%-26s", paste0("p=", PS)), collapse = "")))
  for (k in RULES) {
    v <- vapply(PS, function(pp) {
      w <- S[S$p == pp & S$model == m & S$rule == k, ]
      if (!nrow(w)) "  -" else
        sprintf("%s/%s %s/%s", fm(w$sure4), fm(w$sure20),
                fm(w$ermax, 1), fm(w$medmax, 0))
    }, character(1))
    cat(sprintf("%-20s %s\n", k, paste(sprintf("%-26s", v), collapse = "")))
  }
}
cat("\nformat : Sure-4/Sure-20  E(Rmax)/Med(Rmax)\n")
