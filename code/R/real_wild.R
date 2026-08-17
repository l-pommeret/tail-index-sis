## Application "perplexite en milieu naturel" : les sept regles en concurrence.
##
## Meme protocole que l'application Crime, sur un jeu ou la reponse a une
## justification EVT structurelle : si la surprisal par jeton a une queue
## superieure exponentielle de taux lambda, alors Y = exp(max_t S_t) est
## Pareto d'indice gamma = 1/lambda. Le plateau de Hill le confirme
## (gamma ~ 1.24-1.44 sur k/n dans [0.02, 0.25], QQ exponentiel r = 0.994).
##
## X = profil de surprisal d'un PETIT modele (pythia-70m, 160m) sur le meme
## span de 512 jetons ; Y = exp(max) sous pythia-410m. Aucun quantile de X
## n'est fonction du vecteur de NLL qui definit Y.
##
## Il n'y a pas de verite, mais il y a un CRITERE EXTERNE : la dispersion
## inter-strates de Hill (code/py/gamma_varies_ppl.py), calibree par strates
## aleatoires, qui n'utilise aucun des ecrans compares. On rapporte donc, pour
## chaque regle, son accord de rang avec ce critere -- et le sort qu'elle
## reserve a deux temoins de LOCALISATION PURE (la surprisal moyenne), dont le
## profil de Hill local est plat alors qu'ils deplacent le quantile 0.95.
##
## usage: Rscript code/R/real_wild.R INDIR GAMMAVAR.csv OUT.rds [CORES]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/rank_rules.R")
source("code/R/yoshida_umezu.R")
source("code/R/qa_sis.R")

args <- commandArgs(trailingOnly = TRUE)
INDIR <- if (length(args) >= 1L) args[1L] else "results/wild/ppl100k_ranks"
GVAR <- if (length(args) >= 2L) args[2L] else
  "results/wild/ppl100k/ppl_cov_gammavar.csv"
OUT <- if (length(args) >= 3L) args[3L] else "results/wild/real_wild.rds"
CORES <- if (length(args) >= 4L) as.integer(args[4L]) else 40L

EPS <- 0.05
ASTAR <- 0.35; BSTAR <- 0.15
AGRID <- c(0.30, 0.35, 0.40); BGRID <- c(0.10, 0.15, 0.20)
SEED_TIE <- 973000019L
GAMMA_CARRIERS <- c("s160m_max", "s70m_max", "s160m_kurt", "s70m_kurt")
LOCATION_ONLY <- c("s160m_mean", "s70m_mean")

meta <- readLines(file.path(INDIR, "wild_meta.txt"))
n <- as.integer(meta[1L]); p <- as.integer(meta[2L]); cols <- meta[3:(2L + p)]
Xr <- matrix(readBin(file.path(INDIR, "wild_ranks.bin"), "double", n * p), n, p)
## garde-fou d'ordre memoire (numpy ecrit en C, matrix() lit en Fortran)
for (j in sample.int(p, min(5L, p)))
  stopifnot(length(unique(Xr[, j])) == n, min(Xr[, j]) == 1, max(Xr[, j]) == n)
logY <- readBin(file.path(INDIR, "wild_y.bin"), "double", n); y <- exp(logY)
cat(sprintf("n = %d, p = %d\n", n, p))

## --- les sept regles --------------------------------------------------------
S <- matrix(NA_real_, p, 9L); ii <- 0L
for (a in AGRID) for (b in BGRID) {
  ii <- ii + 1L; al <- n^(-a); hh <- n^(-b) / 2
  S[, ii] <- unlist(parallel::mclapply(seq_len(p), function(j)
    score_coordinate_cpp(Xr[, j], y, hh, al, EPS)[["score"]], mc.cores = CORES))
}
idx <- which(rep(AGRID, each = 3L) == ASTAR & rep(BGRID, 3L) == BSTAR)
u_tb <- tiebreak_u(p, SEED_TIE)
rk <- function(o) { r <- integer(p); r[o] <- seq_len(p); r }
R <- list("tail selected" = rk(order_selected_new(S, idx, u_tb)),
          "tail aggregated" = rk(order_agg_new(S, u_tb)))
cat("  ecran par indice de queue : fait\n"); flush.console()

U <- Xr / (n + 1)
yu <- function(k, h) {
  s <- unlist(parallel::mclapply(seq_len(p), function(j)
    yu_score_coordinate(U[, j], y, k, h)[["score"]], mc.cores = CORES))
  rk(order(-s, u_tb))                       # score grand = deviation locale
}
R[["Yoshida-Umezu paper"]] <- yu(floor(0.072 * n), 1)
R[["Yoshida-Umezu h=2"]] <- yu(floor(0.05 * n), 2)
cat("  Yoshida-Umezu : fait\n"); flush.console()

## qa_sis_scores boucle en serie ; a n = 1e5 chaque rq.fit.br coute plusieurs
## secondes, donc on parallelise sur les coordonnees. Le calcul par coordonnee
## est identique a celui de qa_sis.R -- seule l'iteration change.
qa_par <- function(tau, df = 3L) {
  qm <- as.numeric(quantile(y, tau, type = 7))
  unlist(parallel::mclapply(seq_len(p), function(j) {
    B <- splines::bs(U[, j], df = df)
    fit <- tryCatch(quantreg::rq.fit.br(cbind(1, B), y, tau = tau),
                    error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    mean((cbind(1, B) %*% fit$coefficients - qm)^2)
  }, mc.cores = CORES))
}
for (tau in c(0.90, 0.95, 0.99)) {
  t0 <- proc.time()[3]
  q <- qa_par(tau)
  R[[sprintf("quantile %.2f", tau)]] <- rk(order(-q, u_tb))
  cat(sprintf("    quantile %.2f en %.0f s\n", tau, proc.time()[3] - t0))
  flush.console()
}
cat("  quantile SIS : fait\n"); flush.console()

## --- critere externe --------------------------------------------------------
g <- read.csv(GVAR)
g <- g[match(cols, g$covariate), ]
ok <- !is.na(g$sd_inter)
cat(sprintf("  critere externe disponible pour %d coordonnees sur %d\n",
            sum(ok), p))

## --- tableau de comparaison -------------------------------------------------
rows <- lapply(names(R), function(nm) {
  r <- R[[nm]]
  data.frame(
    rule = nm,
    rho_gamma = round(-cor(r[ok], g$sd_inter[ok], method = "spearman"), 3),
    carriers_top4 = sum(r[match(GAMMA_CARRIERS, cols)] <= 4L),
    carriers_top20 = sum(r[match(GAMMA_CARRIERS, cols)] <= 20L),
    loc_top20 = sum(r[match(LOCATION_ONLY, cols)] <= 20L),
    rank_s160m_max = r[match("s160m_max", cols)],
    rank_s160m_mean = r[match("s160m_mean", cols)],
    rank_s70m_mean = r[match("s70m_mean", cols)])
})
tab <- do.call(rbind, rows)
cat("\n=== comparaison des regles ===\n")
cat("rho_gamma : accord de rang avec le critere externe de variation de gamma\n")
cat("            (positif = la regle classe haut ce que le critere juge variable)\n")
cat("carriers  : les 4 covariables de plus forte variation de gamma\n")
cat("loc_top20 : les 2 temoins de localisation pure dans les 20 premiers\n\n")
print(tab, row.names = FALSE)

cat("\n=== dix premieres par regle ===\n")
for (nm in names(R))
  cat(sprintf("  %-22s %s\n", nm,
              paste(cols[order(R[[nm]])][1:10], collapse = ", ")))

saveRDS(list(n = n, p = p, cols = cols, scores = S, ranks = R,
             gammavar = g, table = tab, astar = ASTAR, bstar = BSTAR),
        OUT, compress = "xz")
write.csv(tab, sub("\\.rds$", "_table.csv", OUT), row.names = FALSE)
cat(sprintf("\nECRIT %s\n", OUT))
