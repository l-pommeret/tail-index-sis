## A partir de quel n chaque regle verrouille-t-elle la cible ?
##
## Le balayage precedent partait de n = 5000, ou toutes les regles reussissent
## deja : il ne discriminait pas. On descend a n = 1000. A ces tailles un seul
## sous-echantillon est trop bruite, donc REP tirages independants par n, et on
## rapporte la moyenne.
##
## La cible reste le consensus des SEPT regles a n = 1e5 (top-10 pour chacune) :
## definie sur le maximum d'information disponible, elle ne privilegie aucune
## regle et ne depend pas des petits echantillons qu'on evalue.
##
## usage: Rscript code/R/wild_smalln.R INDIR TARGET.rds OUTDIR [CORES] [REP]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/rank_rules.R"); source("code/R/yoshida_umezu.R")
suppressPackageStartupMessages({library(quantreg); library(splines)})

args <- commandArgs(trailingOnly = TRUE)
INDIR <- if (length(args) >= 1L) args[1L] else "results/wild/ppl100k_ranks"
TGT <- if (length(args) >= 2L) args[2L] else "results/wild/nsweep_all/nsweep_all.rds"
OUTDIR <- if (length(args) >= 3L) args[3L] else "results/wild/smalln"
CORES <- if (length(args) >= 4L) as.integer(args[4L]) else 40L
REP <- if (length(args) >= 5L) as.integer(args[5L]) else 20L
## Convention de reglage. "fixed_alpha" tient alpha = 0.1 a tous les n : le
## balayage mesure alors l'effet de la taille d'echantillon, a regime constant.
## "fixed_a" tient a = 0.35, la convention publiee : alpha = n^(-0.35) decroit
## avec n, ce qui reproduit ce que ferait un utilisateur appliquant le reglage
## de l'article. Les deux repondent a des questions differentes.
MODE <- if (length(args) >= 6L) args[6L] else "fixed_alpha"
A_FIX <- 0.35
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05
NS <- c(1000L, 1500L, 2000L, 3000L, 5000L, 10000L)

meta <- readLines(file.path(INDIR, "wild_meta.txt"))
N <- as.integer(meta[1L]); p <- as.integer(meta[2L]); cols <- meta[3:(2L + p)]
Xr <- matrix(readBin(file.path(INDIR, "wild_ranks.bin"), "double", N * p), N, p)
for (j in sample.int(p, 5L)) stopifnot(length(unique(Xr[, j])) == N)
logY <- readBin(file.path(INDIR, "wild_y.bin"), "double", N)

a0 <- readRDS(TGT); big <- a0$res[[as.character(max(a0$NS))]]
TARGET <- cols[rowSums(sapply(big, function(r) r <= 10L)) == length(big)]
cat(sprintf("mode = %s ; cible (consensus des %d regles a n = %d) : %s\n",
            MODE, length(big), max(a0$NS), paste(TARGET, collapse = ", ")))
for (nn in NS) { a <- if (MODE == "fixed_a") A_FIX else log(10)/log(nn)
  cat(sprintf("  n = %5d : a = %.3f, alpha = %.4f, 2 n alpha h = %.0f\n",
      nn, a, nn^(-a), 2*nn*nn^(-a)*nn^(-0.15)/2)) }
cat("\n")
tix <- match(TARGET, cols)
u_tb <- tiebreak_u(p, 811000033L)
rk <- function(o) { r <- integer(p); r[o] <- seq_len(p); r }

all_rules <- function(Xs, ys, nn) {
  a <- if (MODE == "fixed_a") A_FIX else log(10) / log(nn)
  AG <- sort(a + c(-0.05, 0, 0.05))
  BG <- c(0.10, 0.15, 0.20)
  S <- vapply(seq_len(9L), function(i) {
    al <- nn^(-AG[(i - 1L) %/% 3L + 1L]); hh <- nn^(-BG[(i - 1L) %% 3L + 1L]) / 2
    unlist(parallel::mclapply(seq_len(p), function(j)
      score_coordinate_cpp(Xs[, j], ys, hh, al, EPS)[["score"]], mc.cores = CORES))
  }, numeric(p))
  U <- Xs / (nn + 1)
  yu <- function(k, h) rk(order(-unlist(parallel::mclapply(seq_len(p), function(j)
    yu_score_coordinate(U[, j], ys, k, h)[["score"]], mc.cores = CORES)), u_tb))
  qa <- function(tau) {
    qm <- as.numeric(quantile(ys, tau, type = 7))
    rk(order(-unlist(parallel::mclapply(seq_len(p), function(j) {
      B <- bs(U[, j], df = 3L)
      f <- tryCatch(rq.fit.br(cbind(1, B), ys, tau = tau), error = function(e) NULL)
      if (is.null(f)) return(NA_real_)
      mean((cbind(1, B) %*% f$coefficients - qm)^2)
    }, mc.cores = CORES)), u_tb))
  }
  list("tail selected" = rk(order_selected_new(S, 5L, u_tb)),
       "tail aggregated" = rk(order_agg_new(S, u_tb)),
       "Yoshida-Umezu paper" = yu(max(2L, floor(0.072 * nn)), 1),
       "Yoshida-Umezu h=2" = yu(max(2L, floor(0.05 * nn)), 2),
       "quantile 0.90" = qa(0.90), "quantile 0.95" = qa(0.95),
       "quantile 0.99" = qa(0.99))
}

rows <- list()
for (nn in NS) for (r in seq_len(REP)) {
  set.seed(90000L + nn + r)
  id <- sample.int(N, nn)
  Xs <- apply(Xr[id, , drop = FALSE], 2L, rank, ties.method = "first")
  R <- all_rules(Xs, exp(logY[id]), nn)
  for (nm in names(R))
    rows[[length(rows) + 1L]] <- data.frame(n = nn, rep = r, regle = nm,
      hits10 = sum(R[[nm]][tix] <= 10L), hits20 = sum(R[[nm]][tix] <= 20L),
      worst = max(R[[nm]][tix]))
  if (r == REP) { cat(sprintf("%s n = %5d fait\n", format(Sys.time(), "%H:%M:%S"), nn)); flush.console() }
}
D <- do.call(rbind, rows)
saveRDS(list(D = D, TARGET = TARGET), file.path(OUTDIR, "smalln.rds"), compress = "xz")
write.csv(D, file.path(OUTDIR, "smalln_raw.csv"), row.names = FALSE)

agg <- function(v) {
  m <- tapply(D[[v]], list(D$regle, D$n), mean)
  round(m[, order(as.integer(colnames(m))), drop = FALSE], 2)
}
cat(sprintf("\n=== cible dans le top-10 (sur %d), moyenne sur %d tirages ===\n",
            length(TARGET), REP)); print(agg("hits10"))
cat("\n=== pire rang parmi la cible, moyenne ===\n"); print(agg("worst"))
write.csv(agg("hits10"), file.path(OUTDIR, "hits10.csv"))
write.csv(agg("worst"), file.path(OUTDIR, "worst.csv"))
