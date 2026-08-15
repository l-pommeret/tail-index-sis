## Application LLM : criblage, calibration par permutation, stabilite.
##
## Reponses (une par colonne du CSV) :
##   Y       = exp(max_t S_t) du petit modele  -- indice de queue de la perte
##   Y_delta = exp(max_t (S_petit - S_grand))  -- exces de perte sur la reference
##
## Regles comparees, comme dans la section de simulation :
##   screen propose au reglage retenu (a*, b*) = (0.30, 0.15)
##   screen agrege par rang minimum sur la grille 3x3
##   quantile SIS a tau = 0.90, 0.95, 0.99
##
## Trois validations que l'application crime ne permettait pas :
##   1. CALIBRATION PAR PERMUTATION. On permute les lignes de X contre Y, ce qui
##      preserve la dependance entre colonnes et detruit le lien a la reponse.
##      B repliques donnent la loi nulle du score minimal, donc un seuil pour
##      "score anormalement bas" -- l'article n'a aujourd'hui aucune calibration
##      inferentielle de Delta_j.
##   2. FREQUENCES DE SELECTION par sous-echantillons disjoints.
##   3. FILTRE DE CONTINUITE strict : au moins 100 valeurs distinctes ET aucun
##      atome portant plus de 2% de la masse. La geometrie des fenetres de rangs
##      en depend, donc la regle de departage des ex aequo est fixee (aleatoire,
##      graine fixee) et non laissee au tri par defaut.
##
## usage: Rscript code/R/screen_llm.R DATA.csv OUTDIR [CORES] [B_PERM]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/qa_sis.R")

args <- commandArgs(trailingOnly = TRUE)
DATA <- args[1L]
OUTDIR <- if (length(args) >= 2L) args[2L] else "results/llm"
CORES <- if (length(args) >= 3L) as.integer(args[3L]) else 92L
BPERM <- if (length(args) >= 4L) as.integer(args[4L]) else 1000L
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
ASTAR <- 0.30; BSTAR <- 0.15
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.10, 0.15, 0.20)
TAUS <- c(0.90, 0.95, 0.99)
SEED_TIES <- 20260815L
NSUB <- 20L                       # sous-echantillons pour les frequences

cat("lecture de", DATA, "\n"); flush.console()
D <- read.csv(DATA, stringsAsFactors = FALSE)
RESP <- strsplit(Sys.getenv("RESP", "max_S_small,max_delta"), ",")[[1]]
RESP <- intersect(RESP, names(D))
meta <- c("doc_id", "url", "dump", "date", "argmax_tok", "top_S",
          "n_tok_doc", "T_scored", "max_S_small", "max_S_large",
          "max_delta", "mean_S_small", "mean_S_large")
cand <- setdiff(names(D), meta)
n <- nrow(D)
cat(sprintf("n = %d documents, %d covariables candidates\n", n, length(cand)))

## ------------------------------------------------- filtre de continuite -----
keep <- vapply(cand, function(v) {
  x <- D[[v]]
  if (!is.numeric(x) || anyNA(x)) return(FALSE)
  if (length(unique(x)) < 100L) return(FALSE)
  max(table(x)) / length(x) <= 0.02
}, logical(1))
cat(sprintf("filtre de continuite : %d retenues sur %d ; ecartees : %s\n",
            sum(keep), length(cand),
            paste(head(cand[!keep], 12), collapse = ", ")))
cand <- cand[keep]
X <- as.matrix(D[, cand, drop = FALSE])

## Departage des ex aequo : aleatoire a graine fixee, applique une fois pour
## toutes. Le score travaille sur des fenetres de rangs, donc un departage
## systematique (ordre d'apparition) creerait une geometrie artificielle.
set.seed(SEED_TIES)
jitter_ranks <- function(x) rank(x, ties.method = "random")
Xr <- apply(X, 2L, jitter_ranks)
storage.mode(Xr) <- "double"

alpha_of <- function(a) n^(-a)
h_of <- function(b) n^(-b) / 2
cat(sprintf("reglage : alpha = %.4f, h = %.4f, 2 n alpha h = %.0f\n",
            alpha_of(ASTAR), h_of(BSTAR),
            2 * n * alpha_of(ASTAR) * h_of(BSTAR)))

score_all <- function(y, a, b) {
  al <- alpha_of(a); hh <- h_of(b)
  unlist(parallel::mclapply(seq_len(ncol(Xr)), function(j)
    score_coordinate_cpp(Xr[, j], y, hh, al, EPS)[["score"]],
    mc.cores = CORES))
}

res <- list()
for (rp in RESP) {
  y <- exp(D[[rp]])                       # Y = exp(max_t S)
  cat("\n=== reponse", rp, "===\n"); flush.console()
  t0 <- proc.time()[3]
  S <- vapply(seq_len(9L), function(i) {
    a <- AGRID[(i - 1L) %/% 3L + 1L]; b <- BGRID[(i - 1L) %% 3L + 1L]
    score_all(y, a, b)
  }, numeric(ncol(Xr)))
  colnames(S) <- sprintf("a%.2f_b%.2f", rep(AGRID, each = 3), rep(BGRID, 3))
  base <- S[, sprintf("a%.2f_b%.2f", ASTAR, BSTAR)]
  R <- apply(S, 2L, rank, ties.method = "first")
  ## Deux agregations. Le rang MINIMUM promeut une coordonnee des qu'un seul
  ## reglage la classe bien : valide derriere une preselection, ou il ne reste
  ## qu'une vingtaine de concurrentes, il promeut les nulles chanceuses quand on
  ## l'applique seul a plusieurs centaines de coordonnees. La MEDIANE est la
  ## regle appropriee ici, et les deux sont rapportees.
  agg <- apply(R, 1L, min)
  agg_med <- apply(R, 1L, median)
  cat(sprintf("  9 passes du score en %.0f s\n", proc.time()[3] - t0))

  qs <- lapply(TAUS, function(tt) qa_sis_scores(Xr / (n + 1), y, tau = tt))
  names(qs) <- paste0("q", TAUS * 100)

  ## --- calibration par permutation ------------------------------------------
  t0 <- proc.time()[3]
  ## Deux lois nulles : celle du MINIMUM sur les p coordonnees, qui donne un
  ## seuil familial, et celle des scores individuels mis en commun, qui donne un
  ## seuil par covariable. La premiere controle le risque de declarer un signal
  ## quand il n'y en a aucun ; la seconde situe une covariable donnee.
  perm <- parallel::mclapply(seq_len(BPERM), function(b) {
    set.seed(10000L + b)
    yp <- y[sample.int(n)]
    vapply(seq_len(ncol(Xr)), function(j)
      score_coordinate_cpp(Xr[, j], yp, h_of(BSTAR), alpha_of(ASTAR),
                           EPS)[["score"]], numeric(1))
  }, mc.cores = CORES)
  null_min <- vapply(perm, min, numeric(1))
  null_pool <- unlist(perm)
  cat(sprintf("  %d permutations en %.0f s\n", BPERM, proc.time()[3] - t0))
  thr05 <- quantile(null_min, 0.05)      # familial
  thr01 <- quantile(null_min, 0.01)
  mrg05 <- quantile(null_pool, 0.05)     # par covariable
  mrg01 <- quantile(null_pool, 0.01)

  ## --- frequences de selection par sous-echantillons -------------------------
  set.seed(4242L)
  folds <- split(sample.int(n), rep_len(seq_len(NSUB), n))
  freq <- Reduce(`+`, parallel::mclapply(folds, function(id) {
    s <- vapply(seq_len(ncol(Xr)), function(j)
      score_coordinate_cpp(Xr[id, j], y[id], h_of(BSTAR),
                           alpha_of(ASTAR), EPS)[["score"]], numeric(1))
    as.integer(rank(s, ties.method = "first") <= 20L)
  }, mc.cores = CORES)) / NSUB

  out <- data.frame(covariate = cand, score = base, rank_base = rank(base),
                    rank_agg_med = rank(agg_med), rank_agg_min = rank(agg),
                    freq_top20 = freq,
                    below_thr05 = base < thr05, below_thr01 = base < thr01,
                    below_mrg05 = base < mrg05, below_mrg01 = base < mrg01,
                    p_marginal = vapply(base, function(v)
                      mean(null_pool <= v), numeric(1)))
  for (k in names(qs))
    out[[paste0("rank_", k)]] <- rank(-qs[[k]], ties.method = "first")
  out <- out[order(out$rank_agg_med), ]
  write.csv(out, file.path(OUTDIR, paste0("screen_", rp, ".csv")), row.names = FALSE)
  saveRDS(list(scores = S, null_min = null_min, null_pool = null_pool,
               thr = c(thr05, thr01, mrg05, mrg01),
               qs = qs, freq = freq, cand = cand),
          file.path(OUTDIR, paste0("screen_", rp, ".rds")), compress = "xz")

  cat(sprintf(paste0("  seuil de permutation : 5%% = %.5f, 1%% = %.5f ;",
                     " score minimal observe = %.5f\n"),
              thr05, thr01, min(base)))
  cat(sprintf(paste0("  seuil par covariable : 5%% = %.5f, 1%% = %.5f\n"),
              mrg05, mrg01))
  cat(sprintf(paste0("  sous le seuil familial : %d a 5%%, %d a 1%% ;",
                     " sous le seuil par covariable : %d a 5%%, %d a 1%%\n"),
              sum(out$below_thr05), sum(out$below_thr01),
              sum(out$below_mrg05), sum(out$below_mrg01)))
  cat("  dix premieres par rang median agrege :\n")
  print(head(out[, c("covariate", "rank_agg_med", "rank_agg_min", "rank_base",
                     "freq_top20",
                     paste0("rank_q", TAUS * 100))], 10), row.names = FALSE)
  res[[rp]] <- out
}
cat("\nECRIT dans", OUTDIR, "\n")
