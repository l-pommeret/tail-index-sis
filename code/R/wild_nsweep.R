## Le signal est-il reel, ou du bruit que le diagnostic de variation de gamma
## aurait ramasse par hasard ?
##
## Sans verite, un critere externe fonde sur la variation de gamma est
## circulaire : il designe ce qu'un estimateur naif de gamma trouve variable,
## puis on juge les ecrans sur leur accord avec lui. Le test non circulaire est
## la STABILITE EN n. Une covariable active doit voir son rang se stabiliser et
## sa significativite croitre avec n ; du bruit produit des rangs erratiques et
## une significativite qui ne progresse pas.
##
## Les sous-echantillons sont EMBOITES (n = 5000 est inclus dans n = 10000,
## etc.) : la comparaison d'un n a l'autre est appariee, on ne compare pas des
## tirages independants.
##
## Le point de fonctionnement est tenu constant : a = log(10)/log(n) donne
## alpha = 0.1 a tous les n. Sinon alpha s'effondre quand n grandit et l'on
## melangerait l'effet de n avec un changement de regime.
##
## usage: Rscript code/R/wild_nsweep.R INDIR OUTDIR [CORES] [BPERM]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/rank_rules.R")

args <- commandArgs(trailingOnly = TRUE)
INDIR <- if (length(args) >= 1L) args[1L] else "results/wild/ppl100k_ranks"
OUTDIR <- if (length(args) >= 2L) args[2L] else "results/wild/nsweep"
CORES <- if (length(args) >= 3L) as.integer(args[3L]) else 40L
BPERM <- if (length(args) >= 4L) as.integer(args[4L]) else 400L
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

EPS <- 0.05; BSTAR <- 0.15
NS <- c(5000L, 10000L, 25000L, 50000L, 100000L)

meta <- readLines(file.path(INDIR, "wild_meta.txt"))
N <- as.integer(meta[1L]); p <- as.integer(meta[2L]); cols <- meta[3:(2L + p)]
Xr <- matrix(readBin(file.path(INDIR, "wild_ranks.bin"), "double", N * p), N, p)
for (j in sample.int(p, 5L))
  stopifnot(length(unique(Xr[, j])) == N, min(Xr[, j]) == 1, max(Xr[, j]) == N)
logY <- readBin(file.path(INDIR, "wild_y.bin"), "double", N)

## sous-echantillons emboites, ordre fixe une fois
set.seed(31415L); perm_rows <- sample.int(N)

res <- list(); sig <- list()
for (nn in NS) {
  id <- perm_rows[seq_len(nn)]
  Xs <- apply(Xr[id, , drop = FALSE], 2L, rank, ties.method = "first")
  ys <- exp(logY[id])
  a <- log(10) / log(nn)                      # alpha = 0.1 exactement
  alpha <- nn^(-a); h <- nn^(-BSTAR) / 2
  AG <- sort(a + c(-0.05, 0, 0.05))
  S <- vapply(seq_len(9L), function(i) {
    al <- nn^(-AG[(i - 1L) %/% 3L + 1L]); hh <- nn^(-c(0.10, 0.15, 0.20)[(i - 1L) %% 3L + 1L]) / 2
    unlist(parallel::mclapply(seq_len(p), function(j)
      score_coordinate_cpp(Xs[, j], ys, hh, al, EPS)[["score"]], mc.cores = CORES))
  }, numeric(p))
  u <- tiebreak_u(p, 811000033L)
  r_sel <- integer(p); r_sel[order_selected_new(S, 5L, u)] <- seq_len(p)
  r_agg <- integer(p); r_agg[order_agg_new(S, u)] <- seq_len(p)

  nullmin <- vapply(parallel::mclapply(seq_len(BPERM), function(b) {
    set.seed(606000L + b)
    yp <- ys[sample.int(nn)]
    vapply(seq_len(p), function(j)
      score_coordinate_cpp(Xs[, j], yp, h, alpha, EPS)[["score"]], numeric(1))
  }, mc.cores = CORES), min, numeric(1))
  thr <- unname(quantile(nullmin, .05))
  nsig <- sum(S[, 5L] < thr)
  res[[as.character(nn)]] <- data.frame(covariate = cols, n = nn,
                                        rank_sel = r_sel, rank_agg = r_agg,
                                        score = S[, 5L], below = S[, 5L] < thr)
  sig[[as.character(nn)]] <- data.frame(n = nn, alpha = alpha,
    two_nah = 2 * nn * alpha * h, thr = thr, min_score = min(S[, 5L]),
    margin = thr - min(S[, 5L]), n_below = nsig)
  cat(sprintf("%s n = %6d : 2nah = %5.0f, seuil %.5f, min %.5f, marge %+.5f, "
              , format(Sys.time(), "%H:%M:%S"), nn, 2*nn*alpha*h, thr,
              min(S[, 5L]), thr - min(S[, 5L])))
  cat(sprintf("%d sous le seuil\n", nsig)); flush.console()
}
R <- do.call(rbind, res); SG <- do.call(rbind, sig)
saveRDS(list(ranks = R, sig = SG, cols = cols), file.path(OUTDIR, "nsweep.rds"),
        compress = "xz")
write.csv(R, file.path(OUTDIR, "nsweep_ranks.csv"), row.names = FALSE)
write.csv(SG, file.path(OUTDIR, "nsweep_sig.csv"), row.names = FALSE)
cat("\n=== significativite par n ===\n"); print(SG, row.names = FALSE, digits = 4)
cat("\nECRIT", OUTDIR, "\n")
