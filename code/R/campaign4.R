## Campagne sur la suite de modèles figée (code/R/generate4.R).
##
##   n   = 2000
##   p   in {500, 1000, 2000}
##   rho in {0, 0.25}
##   modèles M1, M2, M3, M4 (M2 en amas corrélé, M4 à beta = 0.80)
##
## Règles comparées, toutes sur les mêmes jeux de données :
##   1 screen proposé au réglage publié (a*, b*) = (0.30, 0.10)
##   2 screen agrégé par rang minimum sur la grille 3x3
##   3 Yoshida--Umezu, réglage publié h = 1, k = floor(0.072 n)
##   4-7 quantile SIS a tau = 0.90, 0.95, 0.975, 0.99
##   8 pipeline : quantile SIS (tau = 0.95) top-25 puis rang minimum 3x3
##
## Flux de graines 121000007 + cellule*100003 + modele*10007 + r*307, disjoint
## de tous les autres. Cellules ordonnées par coût décroissant, chacune lancée
## comme une seule vague de modèle x réplication sur tous les workers, et
## sauvegardée dès qu'elle est finie (reprise automatique).
##
## usage: Rscript code/R/campaign4.R OUTDIR [NREP] [CORES]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate4.R"); source("code/R/qa_sis.R")
source("code/R/yoshida_umezu.R")

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[1L] else "results/campaign4"
NREP  <- if (length(args) >= 2L) as.integer(args[2L]) else 100L
CORES <- if (length(args) >= 3L) as.integer(args[3L]) else 92L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
ASTAR <- as.numeric(Sys.getenv("ASTAR", "0.30"))
BSTAR <- as.numeric(Sys.getenv("BSTAR", "0.10"))
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.05, 0.10, 0.15)
N <- 2000L
PS <- c(500L, 1000L, 2000L)
RHOS <- c(0, 0.25)
MODELS <- c("M1", "M2", "M3", "M4")
D1 <- 25L
DS <- c(4, 10, 20, 30, 50)
TAUS <- c(0.90, 0.95, 0.975, 0.99)
RULES <- c("screen seul", "screen 9 minimum", "Yoshida-Umezu",
           paste0("quantile .", c("90", "95", "975", "99")),
           "pipeline + 9 min")

cells <- expand.grid(p = PS, rho = RHOS, KEEP.OUT.ATTRS = FALSE)
cells$id <- sprintf("c4_n%04d_p%04d_r%03.0f", N, cells$p, cells$rho * 100)
cells$cost <- 20 * cells$p / 1000 * NREP * length(MODELS)
cells <- cells[order(-cells$cost), ]

score_vec <- function(z, y, p, a, b) {
  h <- N^(-b) / 2; alpha <- N^(-a)
  vapply(seq_len(p), function(j)
    score_coordinate_cpp(z[, j], y, h, alpha, EPS), numeric(5))["score", ]
}

run_one <- function(p, rho, model, r, cellseed) {
  seed <- cellseed + match(model, MODELS) * 10007L + r * 307L
  d <- simulate_dataset4(N, p, rho, model, seed)
  base <- score_vec(d$z, d$y, p, ASTAR, BSTAR)
  R <- matrix(NA_real_, p, length(AGRID) * length(BGRID)); ii <- 0L
  for (a in AGRID) for (b in BGRID) {
    ii <- ii + 1L
    R[, ii] <- rank(score_vec(d$z, d$y, p, a, b), ties.method = "first",
                    na.last = TRUE)
  }
  amin <- apply(R, 1L, min)
  uh <- apply(d$z, 2L, function(x) rank(x, ties.method = "average") / (N + 1))
  yu <- yu_score_matrix(uh, d$y, k = floor(0.072 * N), h = 1)
  ord_yu <- order(-yu$scores, seq_len(p))
  ord_q <- lapply(TAUS, function(tt)
    order(-qa_sis_scores(uh, d$y, tau = tt), seq_len(p)))
  ## l'étage 1 du pipeline utilise tau = 0.95
  surv <- ord_q[[match(0.95, TAUS)]][seq_len(D1)]
  r_pipe <- if (!all(A_GAMMA %in% surv)) {
    Inf
  } else {
    max(match(A_GAMMA, surv[order(amin[surv], surv)]))
  }
  ords <- c(list(order(base, seq_len(p), na.last = TRUE),
                 order(amin, seq_len(p)), ord_yu), ord_q)
  list(model = model, n = N, p = p, rho = rho, replicate = r, seed = seed,
       rmax = c(vapply(ords, function(o) max(match(A_GAMMA, o)), numeric(1)),
                r_pipe),
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
  cellseed <- 121000007L + i * 100003L
  grid <- expand.grid(model = MODELS, r = seq_len(NREP),
                      stringsAsFactors = FALSE)
  out <- parallel::mclapply(seq_len(nrow(grid)), function(k)
    run_one(cells$p[i], cells$rho[i], grid$model[k], grid$r[k], cellseed),
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
ords_n <- seq_len(length(RULES) - 1L)   # toutes les regles sauf le pipeline
files <- sort(list.files(outdir, pattern = "^c4_.*\\.rds$", full.names = TRUE))
S <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f)
  do.call(rbind, lapply(MODELS, function(m) {
    w <- x[vapply(x, function(z) z$model == m, logical(1))]
    do.call(rbind, lapply(seq_along(RULES), function(k) {
      v <- vapply(w, function(z) z$rmax[k], numeric(1))
      r <- data.frame(n = N, p = w[[1]]$p, rho = w[[1]]$rho, model = m,
                      rule = RULES[k], reps = length(w))
      for (d in DS) r[[paste0("sure", d)]] <- mean(v <= d)
      r$ermax <- mean(v[is.finite(v)])
      r$top4_scale <- if (k > length(ords_n)) NA_real_ else
        mean(vapply(w, function(z) z$t4s[k], numeric(1)))
      r$top24_scale <- if (k > length(ords_n)) NA_real_ else
        mean(vapply(w, function(z) z$t24s[k], numeric(1)))
      r
    }))
  }))
}))
write.csv(S, file.path(outdir, "summary.csv"), row.names = FALSE)

fm <- function(v) formatC(v, format = "f", digits = 3)
cat(sprintf("\n=== Suite figee, n=%d, %d replications (erreur type <= %.3f) ===\n",
            N, NREP, 0.5 / sqrt(NREP)))
for (rr in RHOS) {
  cat(sprintf("\n########## rho = %.2f ##########\n", rr))
  for (m in MODELS) {
    cat(sprintf("\n%s   (Sure-4 / Sure-20)\n", m))
    cat(sprintf("%-20s %s\n", "regle",
                paste(sprintf("%-16s", paste0("p=", PS)), collapse = "")))
    for (k in RULES) {
      v <- vapply(PS, function(pp) {
        w <- S[S$p == pp & S$rho == rr & S$model == m & S$rule == k, ]
        if (!nrow(w)) "  -" else sprintf("%s/%s", fm(w$sure4), fm(w$sure20))
      }, character(1))
      cat(sprintf("%-20s %s\n", k, paste(sprintf("%-16s", v), collapse = "")))
    }
  }
}
cat("\n")
