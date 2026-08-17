## Criblage par indice de queue sur la perplexite en milieu naturel.
##
## Quelles caracteristiques de surface d'un document gouvernent la QUEUE de la
## distribution de surprise d'un modele de langue ? La litterature de curation
## par perplexite (Thrush et al., Wu et al. 2025) ne travaille qu'en moyenne
## conditionnelle ; l'indice de queue est le regime que personne n'y regarde.
##
## Y = exp(max_t NLL_t) sur T = 512 jetons apres 64 d'amorce, donc log Y est
## directement l'exces de surprisal. Les covariables sont calculees sur LA MEME
## fenetre, ce qui neutralise la longueur par construction.
##
## Trois blocs de covariables, structure calquee sur la famille A :
##   surface    interpretable, ce qu'on espere voir remonter
##   trigram    frequences de trigrammes denses, semi-interpretable
##   encoder    coordonnees d'encodeur -- le fond de nuisance correle
##
## Regles de classement : celles de l'audit RERUN_RANKFIX (invariantes par
## permutation des colonnes), pas le departage par indice.
##
## usage: Rscript code/R/screen_wild.R INDIR OUTDIR [CORES] [BPERM]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/rank_rules.R")

args <- commandArgs(trailingOnly = TRUE)
INDIR <- if (length(args) >= 1L) args[1L] else "results/wild/ranks"
OUTDIR <- if (length(args) >= 2L) args[2L] else "results/wild/screen"
CORES <- if (length(args) >= 3L) as.integer(args[3L]) else 30L
BPERM <- if (length(args) >= 4L) as.integer(args[4L]) else 1000L
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
## Le reglage a* = 0.35 a ete choisi a n = 2000, ou il donne alpha = 0.070.
## Comme alpha = n^(-a), garder a fixe fait CHUTER alpha quand n grandit :
## 0.018 a n = 1e5, soit un tout autre point de fonctionnement. La theorie
## demande n*alpha*h*Delta^2 >> log(pn), donc un alpha qui ne s'effondre pas.
## WILD_A permet de tenir alpha ~ 0.1 : a = log(10)/log(n).
ASTAR <- as.numeric(Sys.getenv("WILD_A", "0.35"))
BSTAR <- as.numeric(Sys.getenv("WILD_B", "0.15"))
AGRID <- sort(ASTAR + c(-0.05, 0, 0.05)); BGRID <- c(0.10, 0.15, 0.20)
stopifnot(BSTAR %in% BGRID)
SEED_TIE <- 811000033L

base <- file.path(INDIR, "wild")
meta <- readLines(paste0(base, "_meta.txt"))
n <- as.integer(meta[1L]); p <- as.integer(meta[2L])
cols <- meta[3:(2L + p)]
Xr <- matrix(readBin(paste0(base, "_ranks.bin"), "double", n * p), n, p)
logY <- readBin(paste0(base, "_y.bin"), "double", n)
y <- exp(logY)
## Garde-fou : numpy.tofile ecrit en ordre C, matrix() lit en ordre Fortran.
## Une lecture mal ordonnee donne une matrice brouillee et un criblage qui lit
## du bruit -- silencieusement. On verifie donc que chaque colonne est bien une
## permutation de 1..n avant de scorer quoi que ce soit.
chk <- sample.int(p, min(5L, p))
for (j in chk) stopifnot(length(unique(Xr[, j])) == n,
                         min(Xr[, j]) == 1, max(Xr[, j]) == n)
blocks <- read.csv(file.path(INDIR, "wild_blocks.csv"))
stopifnot(nrow(blocks) == p, identical(blocks$name, cols))

alpha <- n^(-ASTAR); h <- n^(-BSTAR) / 2
cat(sprintf("n = %d, p = %d\n", n, p))
cat(sprintf("  blocs : %s\n", paste(sprintf("%s %d", names(table(blocks$block)),
                                            as.integer(table(blocks$block))),
                                    collapse = ", ")))
cat(sprintf("  alpha = %.5f, h = %.5f, 2 n alpha h = %.0f\n",
            alpha, h, 2 * n * alpha * h))
cat(sprintf("  barre de competition a p = %d : %.2f sigma\n", p, qnorm(1 - 1/p)))
cat(sprintf("  log(pn) = %.2f ; n alpha h Delta^2 a Delta=0.186 : %.2f\n",
            log(p * as.double(n)), n * alpha * h * 0.186^2))
flush.console()

score_pass <- function(a, b) {
  al <- n^(-a); hh <- n^(-b) / 2
  unlist(parallel::mclapply(seq_len(p), function(j)
    score_coordinate_cpp(Xr[, j], y, hh, al, EPS)[["score"]], mc.cores = CORES))
}

t0 <- proc.time()[3]
S <- vapply(seq_len(9L), function(i)
  score_pass(AGRID[(i - 1L) %/% 3L + 1L], BGRID[(i - 1L) %% 3L + 1L]),
  numeric(p))
IDX <- which(rep(AGRID, each = 3L) == ASTAR & rep(BGRID, times = 3L) == BSTAR)
stopifnot(length(IDX) == 1L)
cat(sprintf("  9 passes en %.0f s\n", proc.time()[3] - t0)); flush.console()

u <- tiebreak_u(p, SEED_TIE)
ord_sel <- order_selected_new(S, IDX, u)
ord_agg <- order_agg_new(S, u)
rank_sel <- integer(p); rank_sel[ord_sel] <- seq_len(p)
rank_agg <- integer(p); rank_agg[ord_agg] <- seq_len(p)

## --- calibration par permutation --------------------------------------------
## Sans verite sur les covariables reelles, c'est elle qui decide : un rang
## faible ne vaut rien s'il est atteignable par du bruit.
t0 <- proc.time()[3]
perm <- parallel::mclapply(seq_len(BPERM), function(b) {
  set.seed(505000L + b)
  yp <- y[sample.int(n)]
  vapply(seq_len(p), function(j)
    score_coordinate_cpp(Xr[, j], yp, h, alpha, EPS)[["score"]], numeric(1))
}, mc.cores = CORES)
null_min <- vapply(perm, min, numeric(1))
null_pool <- unlist(perm)
cat(sprintf("  %d permutations en %.0f s\n", BPERM, proc.time()[3] - t0))
thr <- c(fam05 = unname(quantile(null_min, .05)),
         fam01 = unname(quantile(null_min, .01)),
         mrg05 = unname(quantile(null_pool, .05)),
         mrg01 = unname(quantile(null_pool, .01)))
pmarg <- vapply(S[, IDX], function(s) mean(null_pool <= s), numeric(1))

## --- frequences de selection par sous-echantillonnage ------------------------
NSUB <- 20L
set.seed(707L)
folds <- split(sample.int(n), rep_len(seq_len(NSUB), n))
freq <- Reduce(`+`, parallel::mclapply(folds, function(id) {
  s <- vapply(seq_len(p), function(j)
    score_coordinate_cpp(Xr[id, j], y[id], h, alpha, EPS)[["score"]], numeric(1))
  as.integer(rank(s, ties.method = "average") <= 20L)
}, mc.cores = CORES)) / NSUB

## --- concurrent : quantile SIS ----------------------------------------------
source("code/R/qa_sis.R")
q95 <- qa_sis_scores(Xr / (n + 1), y, tau = 0.95)
rank_q95 <- rank(-q95, ties.method = "average")

out <- data.frame(covariate = cols, block = blocks$block, score = S[, IDX],
                  rank_sel = rank_sel, rank_agg = rank_agg,
                  rank_q95 = rank_q95, freq_top20 = freq, p_marginal = pmarg,
                  below_fam05 = S[, IDX] < thr["fam05"],
                  below_mrg05 = S[, IDX] < thr["mrg05"],
                  below_mrg01 = S[, IDX] < thr["mrg01"])
out <- out[order(out$rank_sel), ]
write.csv(out, file.path(OUTDIR, "screen_wild.csv"), row.names = FALSE)
saveRDS(list(S = S, null_min = null_min, thr = thr, freq = freq, cols = cols,
             blocks = blocks, logY = logY),
        file.path(OUTDIR, "screen_wild.rds"), compress = "xz")

cat(sprintf("\n  seuil familial 5%% = %.5f ; score minimal observe = %.5f -> %s\n",
            thr["fam05"], min(S[, IDX]),
            if (min(S[, IDX]) < thr["fam05"]) "SIGNIFICATIF" else "non significatif"))
cat(sprintf("  sous seuil familial : %d ; par covariable : %d a 5%% (attendu %.0f), %d a 1%% (attendu %.0f)\n",
            sum(out$below_fam05), sum(out$below_mrg05), 0.05 * p,
            sum(out$below_mrg01), 0.01 * p))
cat("\n  composition des 20 premiers (regle selectionnee) :\n")
print(table(out$block[out$rank_sel <= 20]))
cat("\n  composition des 20 premiers (quantile .95) :\n")
print(table(out$block[out$rank_q95 <= 20]))
cat("\n  quinze premieres (regle selectionnee) :\n")
print(head(out[, c("covariate", "block", "rank_sel", "rank_agg", "rank_q95",
                   "freq_top20", "p_marginal")], 15), row.names = FALSE)
