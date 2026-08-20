## A quel n chaque regle verrouille-t-elle les covariables actives ?
##
## Le piege a eviter : definir les "actives" par la stabilite de NOTRE ecran
## reviendrait a se donner raison. La cible est donc definie SYMETRIQUEMENT --
## les covariables que les SEPT regles placent dans leur top-10 a n = 1e5, le
## plus grand echantillon disponible. Aucune regle n'y est privilegiee, et la
## question devient : laquelle atteint cette cible au plus petit n ?
##
## Sous-echantillons emboites, point de fonctionnement tenu constant
## (alpha = 0.1 a tous les n). Pas de permutations : on compare des RANGS.
##
## usage: Rscript code/R/wild_nsweep_all.R INDIR OUTDIR [CORES]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/rank_rules.R"); source("code/R/yoshida_umezu.R")
suppressPackageStartupMessages({library(quantreg); library(splines)})

args <- commandArgs(trailingOnly = TRUE)
INDIR <- if (length(args) >= 1L) args[1L] else "results/wild/ppl100k_ranks"
OUTDIR <- if (length(args) >= 2L) args[2L] else "results/wild/nsweep_all"
CORES <- if (length(args) >= 3L) as.integer(args[3L]) else 40L
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05; BSTAR <- 0.15
NS <- c(5000L, 10000L, 25000L, 50000L, 100000L)

meta <- readLines(file.path(INDIR, "wild_meta.txt"))
N <- as.integer(meta[1L]); p <- as.integer(meta[2L]); cols <- meta[3:(2L + p)]
Xr <- matrix(readBin(file.path(INDIR, "wild_ranks.bin"), "double", N * p), N, p)
for (j in sample.int(p, 5L)) stopifnot(length(unique(Xr[, j])) == N)
logY <- readBin(file.path(INDIR, "wild_y.bin"), "double", N)
set.seed(31415L); perm_rows <- sample.int(N)
u_tb <- tiebreak_u(p, 811000033L)
rk <- function(o) { r <- integer(p); r[o] <- seq_len(p); r }

all_rules <- function(Xs, ys, nn) {
  a <- log(10) / log(nn); AG <- sort(a + c(-0.05, 0, 0.05))
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
       "Yoshida-Umezu paper" = yu(floor(0.072 * nn), 1),
       "Yoshida-Umezu h=2" = yu(floor(0.05 * nn), 2),
       "quantile 0.90" = qa(0.90), "quantile 0.95" = qa(0.95),
       "quantile 0.99" = qa(0.99))
}

RES <- list()
for (nn in NS) {
  t0 <- proc.time()[3]
  id <- perm_rows[seq_len(nn)]
  Xs <- apply(Xr[id, , drop = FALSE], 2L, rank, ties.method = "first")
  RES[[as.character(nn)]] <- all_rules(Xs, exp(logY[id]), nn)
  cat(sprintf("%s n = %6d fait en %.0f s\n", format(Sys.time(), "%H:%M:%S"),
              nn, proc.time()[3] - t0)); flush.console()
}
saveRDS(list(res = RES, cols = cols, NS = NS), file.path(OUTDIR, "nsweep_all.rds"),
        compress = "xz")

## --- cible symetrique : consensus des sept regles a n maximal ---------------
big <- RES[[as.character(max(NS))]]
inTop <- sapply(big, function(r) r <= 10L)
TARGET <- cols[rowSums(inTop) == length(big)]
cat(sprintf("\ncible = top-10 des SEPT regles a n = %d : %d covariables\n  %s\n",
            max(NS), length(TARGET), paste(TARGET, collapse = ", ")))

tab <- do.call(rbind, lapply(names(big), function(nm)
  data.frame(regle = nm, setNames(as.list(vapply(NS, function(nn)
    sum(RES[[as.character(nn)]][[nm]][match(TARGET, cols)] <= 10L), integer(1))),
    paste0("n", NS)))))
cat(sprintf("\n=== cible atteinte (sur %d) dans le top-10, par n ===\n", length(TARGET)))
print(tab, row.names = FALSE)
worstt <- do.call(rbind, lapply(names(big), function(nm)
  data.frame(regle = nm, setNames(as.list(vapply(NS, function(nn)
    max(RES[[as.character(nn)]][[nm]][match(TARGET, cols)]), integer(1))),
    paste0("n", NS)))))
cat("\n=== pire rang parmi la cible, par n (petit = mieux) ===\n")
print(worstt, row.names = FALSE)
write.csv(tab, file.path(OUTDIR, "target_hits.csv"), row.names = FALSE)
write.csv(worstt, file.path(OUTDIR, "target_worst.csv"), row.names = FALSE)
