## Tables LaTeX et figure de réglage pour la suite de modèles figée.
## Produit manuscript/tables/comparison.tex, manuscript/tables/composition.tex
## et manuscript/figures/tuning_heatmaps.pdf depuis results/campaign5 et
## results/tuning4. Aucun chiffre n'est recopié à la main.
## usage: Rscript code/R/latex4.R

MODELS <- c("M1", "M2", "M3", "M4")
RULES <- c("screen seul", "screen 9 minimum", "Yoshida-Umezu",
           "quantile .90", "quantile .95", "quantile .99")
LAB <- c("Tail-index SIS", "Tail-index SIS, aggregated", "Yoshida--Umezu",
         "Quantile SIS ($\\tau=.90$)", "Quantile SIS ($\\tau=.95$)",
         "Quantile SIS ($\\tau=.99$)")
P0 <- 1000L
fm <- function(v, k = 3) formatC(v, format = "f", digits = k)

files <- sort(list.files("results/campaign5", pattern = "^c5_.*\\.rds$",
                         full.names = TRUE))
cells <- lapply(files, readRDS)
cell <- function(p) {
  for (x in cells) if (x[[1]]$p == p) return(x)
  stop("cellule p=", p, " absente")
}
NREP <- sum(vapply(cell(P0), function(z) z$model == "M1", logical(1)))
SE <- 0.5 / sqrt(NREP)

## ------------------------------------------------------- table comparison --
rows <- character(0)
for (m in MODELS) {
  w <- cell(P0)[vapply(cell(P0), function(z) z$model == m, logical(1))]
  for (k in seq_along(RULES)) {
    v <- vapply(w, function(z) z$rmax[k], numeric(1))
    rows <- c(rows, sprintf("%s & %s & %s & %s & %s & %s\\\\",
      if (k == 1L) m else "", LAB[k], fm(mean(v <= 4)), fm(mean(v <= 20)),
      fm(mean(v), 1), fm(median(v), 0)))
  }
  if (m != "M4") rows <- c(rows, "\\addlinespace")
}
tab <- c("\\begin{table}[t]", "\\centering",
sprintf("\\caption{Method comparison at $n=2000$, $p=%d$, $\\rho=0.25$, on the same %d replications per model. The proposed screen uses the tuning $(a^\\star,b^\\star)=(0.30,0.15)$ of Section~\\ref{sec:tuning}; its aggregated form ranks by the smallest rank attained over the nine tunings $a\\in\\{.30,.35,.40\\}\\times b\\in\\{.10,.15,.20\\}$. Yoshida--Umezu uses its published tuning $h=1$, $k=\\lfloor0.072n\\rfloor$, and the quantile screens cubic $B$-splines with three degrees of freedom. Sure-$d$ is the probability that all four tail-index-active variables rank among the $d$ smallest scores ($d=4$ is exact recovery). Standard errors of the probabilities are at most $%.3f$.}",
        P0, NREP, SE),
"\\label{tab:comparison}", "\\begin{tabular}{llrrrr}", "\\toprule",
"Model & Method & Sure-4 & Sure-20 & $\\mathbb{E}(R_{\\max})$ & $\\mathrm{Med}(R_{\\max})$\\\\",
"\\midrule", rows, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tab, "manuscript/tables/comparison.tex")

## ------------------------------------------------------ table composition --
rows <- character(0)
for (m in MODELS) {
  w <- cell(P0)[vapply(cell(P0), function(z) z$model == m, logical(1))]
  for (k in seq_along(RULES)) {
    rows <- c(rows, sprintf("%s & %s & %s & %s\\\\",
      if (k == 1L) m else "", LAB[k],
      fm(mean(vapply(w, function(z) z$t4g[k], numeric(1))), 2),
      fm(mean(vapply(w, function(z) z$t4s[k], numeric(1))), 2)))
  }
  if (m != "M4") rows <- c(rows, "\\addlinespace")
}
tab <- c("\\begin{table}[t]", "\\centering",
sprintf("\\caption{Composition of the four leading positions, same %d replications as Table~\\ref{tab:comparison}. $\\A_\\gamma=\\{1,2,3,4\\}$ are the tail-index-active coordinates; $\\A_{\\rm scale}=\\{5,\\ldots,24\\}$ are the twenty proxies of the latent scale factor, present in M2 only. Entries are averages of $|\\widehat\\A_4\\cap\\A_\\gamma|$ and $|\\widehat\\A_4\\cap\\A_{\\rm scale}|$ over replications.}",
        NREP),
"\\label{tab:composition}", "\\begin{tabular}{llrr}", "\\toprule",
"Model & Method & $|\\widehat\\A_4\\cap\\A_\\gamma|$ & $|\\widehat\\A_4\\cap\\A_{\\rm scale}|$\\\\",
"\\midrule", rows, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tab, "manuscript/tables/composition.tex")

## ------------------------------------------------------- figure de réglage --
tn <- readRDS("results/tuning4/raw.rds")
G <- tn$grid
S20 <- vapply(MODELS, function(m) {
  w <- tn$jobs[vapply(tn$jobs, function(z) z$model == m, logical(1))]
  rowMeans(vapply(w, function(z) z$rmax, numeric(nrow(G))) <= 20)
}, numeric(nrow(G)))
aa <- sort(unique(G$a)); bb <- sort(unique(G$b))
pdf("manuscript/figures/tuning_heatmaps.pdf", width = 8.6, height = 6.6)
par(mfrow = c(2, 2), mar = c(4.0, 4.0, 2.2, 1.0))
cols <- hcl.colors(20, "YlOrRd", rev = TRUE)
for (i in seq_along(MODELS)) {
  z <- outer(aa, bb, Vectorize(function(x, y)
    S20[G$a == x & G$b == y, i]))
  image(aa, bb, z, zlim = c(0, 1), col = cols, main = MODELS[i],
        xlab = "tail exponent a", ylab = "bandwidth exponent b")
  for (u in seq_along(aa)) for (v in seq_along(bb))
    text(aa[u], bb[v], sprintf("%.2f", z[u, v]), cex = .6)
}
dev.off()

## ---------------------------------------------------- figure sensibilite tau --
TAUS <- c(0.90, 0.95, 0.99)
K_TAU <- 4:6; K_AGG <- 2L
pdf("manuscript/figures/tau_sensitivity.pdf", width = 5.8, height = 4.4)
par(mar = c(4.2, 4.4, 1.0, 1.0))
plot(range(TAUS), c(0, 1), type = "n", xlab = expression(tau),
     ylab = "Sure-20 of quantile SIS")
pchs <- c(M1 = 19, M2 = 17, M3 = 15, M4 = 18)
for (m in MODELS) {
  w <- cell(P0)[vapply(cell(P0), function(z) z$model == m, logical(1))]
  v <- vapply(K_TAU, function(k)
    mean(vapply(w, function(z) z$rmax[k], numeric(1)) <= 20), numeric(1))
  lines(TAUS, v, type = "b", pch = pchs[m])
  agg <- mean(vapply(w, function(z) z$rmax[K_AGG], numeric(1)) <= 20)
  points(0.9955, agg, pch = pchs[m], col = "grey55", xpd = NA)
}
legend("bottomleft", legend = MODELS, pch = pchs, lty = 1, bty = "n")
dev.off()

cat("ecrit : comparison.tex, composition.tex, tuning_heatmaps.pdf, tau_sensitivity.pdf\n")
cat(sprintf("controle 2 n alpha h a (a*,b*)=(0.30,0.15) : %.1f\n",
            2 * 2000 * 2000^(-0.30) * 2000^(-0.15) / 2))
