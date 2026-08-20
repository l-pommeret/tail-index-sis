## Ou, dans le plan (n, a), l'ecran par indice de queue domine-t-il ?
##
## Les cinq concurrents ne dependent pas de a : on les calcule UNE fois par
## tirage, et on balaie a pour le seul ecran de queue. La comparaison est donc
## appariee -- meme sous-echantillon, meme cible, seul le reglage change.
##
## Cible : consensus des sept regles a n = 1e5 (top-10 pour chacune), definie
## sur le maximum d'information et ne privilegiant aucune regle.
##
## Critere : "pire rang parmi la cible", petit = mieux. L'avantage rapporte est
##   avantage = (meilleur pire-rang des concurrents) - (pire-rang de la queue)
## positif = la queue fait mieux que le meilleur de ses cinq concurrents.
##
## usage: Rscript code/R/wild_grid_na.R INDIR TARGET.rds OUTDIR [CORES]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/rank_rules.R"); source("code/R/yoshida_umezu.R")
suppressPackageStartupMessages({library(quantreg); library(splines)})

args <- commandArgs(trailingOnly = TRUE)
INDIR <- if (length(args) >= 1L) args[1L] else "results/wild/ppl100k_ranks"
TGT <- if (length(args) >= 2L) args[2L] else "results/wild/nsweep_all/nsweep_all.rds"
OUTDIR <- if (length(args) >= 3L) args[3L] else "results/wild/grid_na"
CORES <- if (length(args) >= 4L) as.integer(args[4L]) else 40L
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
NS <- c(1000L, 2000L, 5000L, 10000L, 25000L, 50000L)
REPS <- c(20L, 20L, 20L, 10L, 5L, 3L)          # le cout est domine par quantile SIS
AS <- c(0.15, 0.20, 0.25, 0.30, 0.35, 0.40)

meta <- readLines(file.path(INDIR, "wild_meta.txt"))
N <- as.integer(meta[1L]); p <- as.integer(meta[2L]); cols <- meta[3:(2L + p)]
Xr <- matrix(readBin(file.path(INDIR, "wild_ranks.bin"), "double", N * p), N, p)
for (j in sample.int(p, 5L)) stopifnot(length(unique(Xr[, j])) == N)
logY <- readBin(file.path(INDIR, "wild_y.bin"), "double", N)
a0 <- readRDS(TGT); big <- a0$res[[as.character(max(a0$NS))]]
TARGET <- cols[rowSums(sapply(big, function(r) r <= 10L)) == length(big)]
tix <- match(TARGET, cols)
cat(sprintf("cible : %s\n\n", paste(TARGET, collapse = ", ")))
u_tb <- tiebreak_u(p, 811000033L)
rk <- function(o) { r <- integer(p); r[o] <- seq_len(p); r }

rows <- list()
for (k in seq_along(NS)) {
  nn <- NS[k]
  for (r in seq_len(REPS[k])) {
    set.seed(120000L + nn + r)
    id <- sample.int(N, nn)
    Xs <- apply(Xr[id, , drop = FALSE], 2L, rank, ties.method = "first")
    ys <- exp(logY[id]); U <- Xs / (nn + 1)

    ## --- concurrents, une seule fois -------------------------------------
    yu <- function(kk, h) rk(order(-unlist(parallel::mclapply(seq_len(p), function(j)
      yu_score_coordinate(U[, j], ys, kk, h)[["score"]], mc.cores = CORES)), u_tb))
    qa <- function(tau) {
      qm <- as.numeric(quantile(ys, tau, type = 7))
      rk(order(-unlist(parallel::mclapply(seq_len(p), function(j) {
        B <- bs(U[, j], df = 3L)
        f <- tryCatch(rq.fit.br(cbind(1, B), ys, tau = tau), error = function(e) NULL)
        if (is.null(f)) return(NA_real_)
        mean((cbind(1, B) %*% f$coefficients - qm)^2)
      }, mc.cores = CORES)), u_tb))
    }
    comp <- list(yu(max(2L, floor(0.072 * nn)), 1), yu(max(2L, floor(0.05 * nn)), 2),
                 qa(0.90), qa(0.95), qa(0.99))
    comp_worst <- min(vapply(comp, function(z) max(z[tix]), numeric(1)))

    ## --- l'ecran de queue, pour chaque a ---------------------------------
    for (a in AS) {
      AG <- sort(a + c(-0.05, 0, 0.05)); BG <- c(0.10, 0.15, 0.20)
      S <- vapply(seq_len(9L), function(i) {
        al <- nn^(-AG[(i - 1L) %/% 3L + 1L]); hh <- nn^(-BG[(i - 1L) %% 3L + 1L]) / 2
        unlist(parallel::mclapply(seq_len(p), function(j)
          score_coordinate_cpp(Xs[, j], ys, hh, al, EPS)[["score"]], mc.cores = CORES))
      }, numeric(p))
      r_sel <- rk(order_selected_new(S, 5L, u_tb))
      r_agg <- rk(order_agg_new(S, u_tb))
      rows[[length(rows) + 1L]] <- data.frame(
        n = nn, rep = r, a = a, alpha = nn^(-a),
        two_nah = 2 * nn * nn^(-a) * nn^(-0.15) / 2,
        agg_worst = max(r_agg[tix]), sel_worst = max(r_sel[tix]),
        agg_hits = sum(r_agg[tix] <= 10L), sel_hits = sum(r_sel[tix] <= 10L),
        comp_worst = comp_worst)
    }
  }
  cat(sprintf("%s n = %6d (%d tirages) fait\n", format(Sys.time(), "%H:%M:%S"),
              nn, REPS[k])); flush.console()
}
D <- do.call(rbind, rows)
saveRDS(list(D = D, TARGET = TARGET), file.path(OUTDIR, "grid_na.rds"), compress = "xz")
write.csv(D, file.path(OUTDIR, "grid_na_raw.csv"), row.names = FALSE)

mat <- function(v) {
  m <- tapply(D[[v]], list(D$a, D$n), mean)
  round(m[, order(as.integer(colnames(m))), drop = FALSE], 2)
}
cat("\n=== pire rang, regle AGREGEE (lignes = a, colonnes = n) ===\n"); print(mat("agg_worst"))
cat("\n=== pire rang, meilleur des CINQ concurrents ===\n")
print(round(tapply(D$comp_worst, D$n, mean)[order(as.integer(names(tapply(D$comp_worst, D$n, mean))))], 2))
cat("\n=== AVANTAGE = concurrents - queue agregee (positif = la queue gagne) ===\n")
D$adv <- D$comp_worst - D$agg_worst
print(mat("adv"))
cat("\n=== 2 n alpha h correspondant ===\n"); print(mat("two_nah"))
write.csv(mat("adv"), file.path(OUTDIR, "advantage.csv"))
write.csv(mat("agg_worst"), file.path(OUTDIR, "agg_worst.csv"))
