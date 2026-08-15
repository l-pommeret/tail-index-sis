## Criblage des etats internes : quelles directions du flux residuel annoncent
## un echec catastrophique a venir ?
##
## X = activations moyennees sur le prefixe (passe), Y = exp(max_t S) sur le
## suffixe (futur). Le decoupage evite la degenerescence que le point 5 du
## protocole signale quand X est une fonction deterministe de la passe qui
## produit Y.
##
## La calibration par permutation n'est pas decorative ici : contrairement au
## bloc longueur, dont le statut est etabli analytiquement, on n'a aucune verite
## sur les directions internes. C'est elle qui decide.
##
## usage: Rscript code/R/screen_internal.R INDIR CONFIG OUTDIR [CORES] [BPERM] [QSIS]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)

args <- commandArgs(trailingOnly = TRUE)
INDIR <- args[1L]; CONFIG <- args[2L]
OUTDIR <- if (length(args) >= 3L) args[3L] else "results/llm/screen_internal"
CORES <- if (length(args) >= 4L) as.integer(args[4L]) else 88L
BPERM <- if (length(args) >= 5L) as.integer(args[5L]) else 1000L
QSIS <- if (length(args) >= 6L) as.logical(args[6L]) else FALSE
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
ASTAR <- 0.30; BSTAR <- 0.15
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.10, 0.15, 0.20)

base <- file.path(INDIR, CONFIG)
meta <- readLines(paste0(base, "_meta.txt"))
n <- as.integer(meta[1L]); p <- as.integer(meta[2L])
cols <- meta[3:(2L + p)]
Xr <- matrix(readBin(paste0(base, "_ranks.bin"), "double", n * p), n, p)
y <- exp(readBin(paste0(base, "_y.bin"), "double", n))
cat(sprintf("config %s : n = %d, p = %d\n", CONFIG, n, p))

alpha <- n^(-ASTAR); h <- n^(-BSTAR) / 2
cat(sprintf("alpha = %.5f, h = %.5f, 2 n alpha h = %.0f, n/p = %.1f\n",
            alpha, h, 2 * n * alpha * h, n / p))
cat(sprintf("barre de competition a p = %d : %.2f sigma\n", p, qnorm(1 - 1 / p)))

score_pass <- function(a, b) {
  al <- n^(-a); hh <- n^(-b) / 2
  unlist(parallel::mclapply(seq_len(p), function(j)
    score_coordinate_cpp(Xr[, j], y, hh, al, EPS)[["score"]], mc.cores = CORES))
}

t0 <- proc.time()[3]
S <- vapply(seq_len(9L), function(i)
  score_pass(AGRID[(i - 1L) %/% 3L + 1L], BGRID[(i - 1L) %% 3L + 1L]),
  numeric(p))
## la grille est parcourue a-externe / b-interne : l'indice de (0.30, 0.15) est 2
IDX <- which(rep(AGRID, each = 3L) == ASTAR & rep(BGRID, times = 3L) == BSTAR)
stopifnot(length(IDX) == 1L)
base_sc <- S[, IDX]
R <- apply(S, 2L, rank, ties.method = "first")
agg_med <- apply(R, 1L, median)
cat(sprintf("  9 passes en %.0f s\n", proc.time()[3] - t0)); flush.console()

## --- calibration par permutation -------------------------------------------
t0 <- proc.time()[3]
perm <- parallel::mclapply(seq_len(BPERM), function(b) {
  set.seed(70000L + b)
  yp <- y[sample.int(n)]
  vapply(seq_len(p), function(j)
    score_coordinate_cpp(Xr[, j], yp, h, alpha, EPS)[["score"]], numeric(1))
}, mc.cores = CORES)
null_min <- vapply(perm, min, numeric(1))
null_pool <- unlist(perm)
cat(sprintf("  %d permutations en %.0f s\n", BPERM, proc.time()[3] - t0))
## unname : quantile() renvoie un vecteur nomme, sinon les noms deviennent
## "fam05.5%" et l'indexation par thr["fam05"] donne NA
thr <- c(fam05 = unname(quantile(null_min, .05)),
         fam01 = unname(quantile(null_min, .01)),
         mrg05 = unname(quantile(null_pool, .05)),
         mrg01 = unname(quantile(null_pool, .01)))

## --- frequences de selection ------------------------------------------------
NSUB <- 20L
set.seed(909L)
folds <- split(sample.int(n), rep_len(seq_len(NSUB), n))
freq <- Reduce(`+`, parallel::mclapply(folds, function(id) {
  s <- vapply(seq_len(p), function(j)
    score_coordinate_cpp(Xr[id, j], y[id], n^(-BSTAR) / 2, n^(-ASTAR),
                         EPS)[["score"]], numeric(1))
  as.integer(rank(s, ties.method = "first") <= 20L)
}, mc.cores = CORES)) / NSUB

out <- data.frame(covariate = cols, score = base_sc,
                  rank_base = rank(base_sc), rank_agg_med = rank(agg_med),
                  freq_top20 = freq,
                  below_fam05 = base_sc < thr["fam05"],
                  below_mrg05 = base_sc < thr["mrg05"],
                  below_mrg01 = base_sc < thr["mrg01"])
if (QSIS) {
  source("code/R/qa_sis.R")
  q95 <- qa_sis_scores(Xr / (n + 1), y, tau = 0.95)
  out$rank_q95 <- rank(-q95, ties.method = "first")
}
out <- out[order(out$rank_base), ]
write.csv(out, file.path(OUTDIR, paste0(CONFIG, ".csv")), row.names = FALSE)
saveRDS(list(scores = S, null_min = null_min, thr = thr, freq = freq,
             cols = cols), file.path(OUTDIR, paste0(CONFIG, ".rds")),
        compress = "xz")

cat(sprintf("  seuil familial 5%% = %.5f ; score minimal observe = %.5f -> %s\n",
            thr["fam05"], min(base_sc),
            if (min(base_sc) < thr["fam05"]) "SIGNIFICATIF" else "non significatif"))
cat(sprintf("  sous seuil familial : %d ; sous seuil par covariable : %d a 5%%, %d a 1%%\n",
            sum(out$below_fam05), sum(out$below_mrg05), sum(out$below_mrg01)))
cat(sprintf("  attendu par hasard au seuil par covariable : %.0f a 5%%, %.0f a 1%%\n",
            0.05 * p, 0.01 * p))
cat(sprintf("  frequence de selection : med %.2f, max %.2f, nb a 1.00 : %d\n",
            median(freq), max(freq), sum(freq == 1)))
cat("  dix premieres :\n")
print(head(out, 10), row.names = FALSE)
