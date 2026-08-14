## Phase de réglage sur la suite de modèles figée (code/R/generate4.R).
##
## Sure-20 sur la grille croisée du papier,
##   a in {.25,.30,.35,.40,.45,.50},  b in {0,.05,.10,.15,.20,.30,.40},
## pour les quatre modèles au design de référence n=2000, p=1000, rho=0.25,
## 200 réplications par cellule.
##
## Chaque réplication génère UN jeu de données et le score aux 42 réglages :
## les cellules du heatmap sont donc appariées, ce qui réduit fortement le
## bruit sur les comparaisons entre réglages, celles-là mêmes qui décident du
## choix de (a*, b*).  Le prix est que les 42 cellules ne sont plus
## indépendantes, ce dont il faut tenir compte en lisant les écarts absolus.
##
## Flux de graines 141000003 + modele*100003 + r*307, disjoint des autres.
##
## usage: Rscript code/R/tuning4.R OUT.rds [NREP] [CORES]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate4.R")

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1L) args[1L] else "results/tuning4/raw.rds"
NREP  <- if (length(args) >= 2L) as.integer(args[2L]) else 200L
CORES <- if (length(args) >= 3L) as.integer(args[3L]) else 92L
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
N <- 2000L; P <- 1000L; RHO <- 0.25
MODELS <- c("M1", "M2", "M3", "M4")
AS <- c(.25, .30, .35, .40, .45, .50)
BS <- c(0, .05, .10, .15, .20, .30, .40)
GRID <- expand.grid(b = BS, a = AS)          # a externe, b interne

jobs <- expand.grid(m = seq_along(MODELS), r = seq_len(NREP))

one <- function(ix) {
  m <- MODELS[jobs$m[ix]]; r <- jobs$r[ix]
  seed <- 141000003L + jobs$m[ix] * 100003L + r * 307L
  d <- simulate_dataset4(N, P, RHO, m, seed)
  rmax <- vapply(seq_len(nrow(GRID)), function(g) {
    h <- N^(-GRID$b[g]) / 2; alpha <- N^(-GRID$a[g])
    s <- vapply(seq_len(P), function(j)
      score_coordinate_cpp(d$z[, j], d$y, h, alpha, EPS), numeric(5))["score", ]
    max(match(A_GAMMA, order(s, seq_len(P), na.last = TRUE)))
  }, numeric(1))
  list(model = m, replicate = r, seed = seed, rmax = rmax)
}

res <- parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = CORES,
                          mc.preschedule = FALSE)
bad <- vapply(res, function(z) inherits(z, "try-error") || is.null(z), logical(1))
if (any(bad)) stop(sum(bad), " echecs : ", as.character(res[bad][[1]]))
saveRDS(list(grid = GRID, models = MODELS, n = N, p = P, rho = RHO,
             nrep = NREP, jobs = res), out_path, compress = "xz")

## ------------------------------------------------------------------ tables --
S20 <- vapply(MODELS, function(m) {
  w <- res[vapply(res, function(z) z$model == m, logical(1))]
  R <- vapply(w, function(z) z$rmax, numeric(nrow(GRID)))
  rowMeans(R <= 20)
}, numeric(nrow(GRID)))
out <- cbind(GRID, as.data.frame(S20))
out$mean <- rowMeans(S20)
write.csv(out, sub("\\.rds$", "_summary.csv", out_path), row.names = FALSE)

fm <- function(v) formatC(v, format = "f", digits = 3)
cat(sprintf("\nReglage sur la suite figee : n=%d p=%d rho=%.2f, %d replications\n",
            N, P, RHO, NREP))
cat(sprintf("erreur type <= %.3f ; comparaisons entre reglages appariees\n\n",
            0.5 / sqrt(NREP)))
for (m in c(MODELS, "mean")) {
  cat(sprintf("%s : Sure-20\n", if (m == "mean") "MOYENNE DES 4 MODELES" else m))
  cat(sprintf("%6s %s\n", "a\\b", paste(sprintf("%7s", BS), collapse = "")))
  for (a in AS) {
    v <- vapply(BS, function(b) out[[m]][out$a == a & out$b == b], numeric(1))
    cat(sprintf("%6.2f %s\n", a, paste(sprintf("%7s", fm(v)), collapse = "")))
  }
  cat("\n")
}
best <- out[which.max(out$mean), ]
cat(sprintf("Meilleur reglage sur la moyenne des quatre modeles : a=%.2f b=%.2f (%.3f)\n",
            best$a, best$b, best$mean))
cat("Meilleur par modele :\n")
for (m in MODELS) {
  bb <- out[which.max(out[[m]]), ]
  cat(sprintf("  %s : a=%.2f b=%.2f (%.3f)\n", m, bb$a, bb$b, bb[[m]]))
}
