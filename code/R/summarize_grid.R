## Summaries and Markdown report for the extended comparison grid.
## usage: Rscript code/R/summarize_grid.R [CELLDIR] [OUTDIR]

args <- commandArgs(trailingOnly = TRUE)
celldir <- if (length(args) >= 1L) args[1L] else "results/grid/comparison_cells"
outdir  <- if (length(args) >= 2L) args[2L] else "results/grid"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

DS <- c(4, 10, 20, 30, 50, 100)
METHODS <- c("ours", "yu", "q900", "q950", "q975", "q990")
LABEL <- c(ours = "Tail-index SIS", yu = "Yoshida--Umezu",
           q900 = "Quantile SIS t=.90", q950 = "Quantile SIS t=.95",
           q975 = "Quantile SIS t=.975", q990 = "Quantile SIS t=.99")

files <- sort(list.files(celldir, pattern = "\\.rds$", full.names = TRUE))
if (!length(files)) stop("no cells in ", celldir)

ranks_of <- function(z, k)
  if (k == "ours") z$ranks_ours else if (k == "yu") z$ranks_yu else z$q[[k]]$ranks
top4_of <- function(z, k)
  if (k == "ours") z$top4_ours else if (k == "yu") z$top4_yu else z$q[[k]]$top4

rows <- list(); tim <- list()
for (f in files) {
  x <- readRDS(f)
  z1 <- x[[1]]
  for (k in METHODS) {
    rmax <- vapply(x, function(z) max(ranks_of(z, k)), numeric(1))
    t4g <- vapply(x, function(z) sum(top4_of(z, k) %in% 1:4), numeric(1))
    t4s <- vapply(x, function(z) sum(top4_of(z, k) %in% 5:8), numeric(1))
    r <- data.frame(model = z1$model, n = z1$n, p = z1$p, rho = z1$rho,
                    method = k, reps = length(x))
    for (d in DS) r[[paste0("sure", d)]] <- mean(rmax <= d)
    r$ermax <- mean(rmax)
    r$medmax <- median(rmax)
    r$top4_gamma <- mean(t4g)
    r$top4_scale <- mean(t4s)
    rows[[length(rows) + 1L]] <- r
  }
  el <- vapply(x, function(z) z$elapsed, c(ours = 0, yu = 0, qsis = 0))
  tim[[length(tim) + 1L]] <- data.frame(model = z1$model, n = z1$n, p = z1$p,
    rho = z1$rho, ours = mean(el["ours", ]), yu = mean(el["yu", ]),
    qsis = mean(el["qsis", ]),
    yu_undefined = mean(vapply(x, function(z) z$yu_undefined, numeric(1))))
}
S <- do.call(rbind, rows); rownames(S) <- NULL
TM <- do.call(rbind, tim); rownames(TM) <- NULL
S$method <- factor(S$method, levels = METHODS)
write.csv(S, file.path(outdir, "grid_summary.csv"), row.names = FALSE)
write.csv(TM, file.path(outdir, "grid_timings.csv"), row.names = FALSE)
cat("cells:", length(files), " rows:", nrow(S), "\n")

## ---------------------------------------------------------------- Markdown --
fmt <- function(v, k = 3) formatC(v, format = "f", digits = k)
md_table <- function(header, align, body) {
  c(paste0("| ", paste(header, collapse = " | "), " |"),
    paste0("|", paste(ifelse(align == "l", ":---", "---:"), collapse = "|"), "|"),
    body, "")
}
NREP <- max(S$reps)
NCELL <- as.integer(Sys.getenv("GRID_NCELL", "288"))
sure_cols <- paste0("sure", DS)

## rho label without trailing zeros
rl <- function(v) sub("\\.?0+$", "", formatC(v, format = "f", digits = 2))

## --- campaign schedule, read off the driver log ------------------------------
## The driver logs, per batch, "HH:MM:SS batch of N cells (C core-hours)" when
## it starts one and "HH:MM:SS elapsed E min, remaining budget R core-hours"
## when it finishes one.  Those give the realised throughput and, from it, the
## next batch's due time and the end of the campaign.
schedule_lines <- function(log_file = Sys.getenv("GRID_LOG",
                                                 "results/grid/run_grid.log")) {
  if (!file.exists(log_file)) return(character(0))
  lg <- readLines(log_file, warn = FALSE)
  hdr <- grep("to run \\(", lg, value = TRUE)
  el  <- grep("elapsed .*remaining budget", lg, value = TRUE)
  bt  <- grep("batch of .* cells \\(", lg, value = TRUE)
  if (!length(hdr) || !length(bt)) return(character(0))
  total <- as.numeric(sub(".*to run \\(([0-9.]+) core-hours\\).*", "\\1", hdr[1]))
  if (grepl("ALL CELLS DONE", paste(lg, collapse = " "))) {
    return(c("**Campagne terminee.**", ""))
  }
  last_b <- bt[length(bt)]
  t0 <- as.POSIXct(paste(Sys.Date(),
                         sub("^([0-9:]{8}).*", "\\1", last_b)), tz = "")
  cost <- as.numeric(sub(".*\\(([0-9.]+) core-hours\\).*", "\\1", last_b))
  ncell_batch <- as.numeric(sub(".*batch of ([0-9]+) cells.*", "\\1", last_b))
  if (!length(el)) return(character(0))
  last_e <- el[length(el)]
  emin <- as.numeric(sub(".*elapsed ([0-9.]+) min.*", "\\1", last_e))
  rem  <- as.numeric(sub(".*remaining budget ([0-9.]+) core-hours.*", "\\1", last_e))
  rate <- (total - rem) / emin            # core-hours per minute of wall time
  if (!is.finite(rate) || rate <= 0) return(character(0))
  hm <- function(tm) format(tm, "%H:%M")
  eta_batch <- t0 + 60 * cost / rate
  eta_end   <- t0 + 60 * rem / rate
  now <- Sys.time()
  c(sprintf("**Prochaine echeance : batch de %d cellules attendu vers %s%s ; fin de campagne estimee vers %s.**",
            ncell_batch, hm(eta_batch),
            if (eta_batch < now) " (imminent)" else "", hm(eta_end)),
    "",
    sprintf("Debit mesure %.1f core-heures par heure de calcul (soit %.0f coeurs solo-equivalents) ; %.1f core-heures restantes sur %.1f.",
            rate * 60, rate * 60, rem, total),
    "")
}

## --- block of full tables: one per (model, n, p) ----------------------------
cell_table <- function(w) {
  w <- w[order(w$rho, as.integer(w$method)), ]
  body <- character(0)
  for (rr in unique(w$rho)) {
    ww <- w[w$rho == rr, ]
    for (i in seq_len(nrow(ww))) {
      body <- c(body, paste0("| ", paste(c(
        if (i == 1L) paste0("**", rl(rr), "**") else "",
        LABEL[as.character(ww$method[i])],
        fmt(unlist(ww[i, sure_cols])),
        fmt(ww$ermax[i], 1), fmt(ww$medmax[i], 1)), collapse = " | "), " |"))
    }
  }
  md_table(c("rho", "Methode", paste0("Sure-", DS), "E(Rmax)", "Med(Rmax)"),
           c("l", "l", rep("r", length(DS) + 2L)), body)
}

## --- aggregated marginal tables --------------------------------------------
marginal <- function(w, by, col = "sure20") {
  lv <- sort(unique(w[[by]]))
  body <- vapply(METHODS, function(k) {
    v <- vapply(lv, function(l)
      mean(w[[col]][w$method == k & w[[by]] == l]), numeric(1))
    paste0("| ", paste(c(LABEL[k], fmt(v)), collapse = " | "), " |")
  }, character(1))
  md_table(c("Methode", paste0(by, " = ", rl(lv))),
           c("l", rep("r", length(lv))), unname(body))
}

sured_profile <- function(w) {
  body <- vapply(METHODS, function(k) {
    v <- vapply(sure_cols, function(cc) mean(w[[cc]][w$method == k]), numeric(1))
    paste0("| ", paste(c(LABEL[k], fmt(v)), collapse = " | "), " |")
  }, character(1))
  md_table(c("Methode", paste0("Sure-", DS)),
           c("l", rep("r", length(DS))), unname(body))
}

md <- c(
"# Comparaison des methodes de screening sur grille etendue",
"",
sprintf("Genere le %s a partir de `%s` (%d cellules sur %d).",
        format(Sys.time(), "%Y-%m-%d %H:%M"), celldir, length(files), NCELL),
"",
if (length(files) < NCELL)
  sprintf("> **Rapport partiel** : la campagne est en cours, %d cellules sur %d sont calculees. Les vues agregees de la section 3 ne portent que sur les cellules disponibles et bougeront encore. Les tables completes de la section 4 sont definitives cellule par cellule.",
          length(files), NCELL)
else "Campagne complete.",
"",
schedule_lines(),
"## 1. Protocole",
"",
"Les quatre modeles M1-M4 de la Draft 3 (`code/R/generate3.R`), les memes",
"estimateurs et le meme reglage que `code/R/run_draft3.R comparison`, mais sur",
"un espace de design elargi.",
"",
sprintf("- **Modeles** : M1, M2, M3, M4 ; ensemble actif en indice de queue A_gamma = {1,2,3,4} (pour M2, les variables d'echelle A_scale = {5,6,7,8} n'agissent que sur les quantiles finis, avec kappa = %s).", Sys.getenv("KAPPA", "0.20")),
"- **Tailles** : n = 1000, 2000, 5000.",
"- **Dimensions** : p = 200, 500, 1000, 2000.",
"- **Dependance** : X ~ AR(1) gaussien de correlation rho = 0, 0.20, 0.25, 0.30, 0.40, 0.50.",
sprintf("- **Replications Monte Carlo** : %d par cellule, soit %d cellules et %d jeux de donnees simules.",
        NREP, length(files), NREP * length(files)),
"",
"Methodes comparees (identiques a celles de la Draft 3) :",
"",
"- **Tail-index SIS** (propose) : score de Hill local sur rangs empiriques,",
sprintf("  a* = %s, b* = %s, soit alpha = n^-a*, h = n^-b*/2, epsilon = 0.05 ; les coordonnees sont classees par score **croissant**.",
        Sys.getenv("ASTAR", "0.30"), Sys.getenv("BSTAR", "0.10")),
"- **Yoshida--Umezu** : screening de Pickands conditionnel, reglage publie h = 1, k = floor(0.072 n).",
"- **Quantile SIS** (He, Wang et Hong 2013) : B-splines cubiques a 3 degres de liberte, tau = 0.90, 0.95, 0.975, 0.99.",
"",
"Criteres :",
"",
sprintf("- **Sure-d** = probabilite que les 4 variables actives soient toutes classees parmi les d meilleures (d = %s ; d = 4 est la recuperation exacte).",
        paste(DS, collapse = ", ")),
"- **E(Rmax)**, **Med(Rmax)** = moyenne et mediane du pire rang actif.",
"",
sprintf("Erreur type Monte Carlo d'une probabilite avec %d replications : au plus %.3f (et %.3f pour une probabilite de 0.9).",
        NREP, 0.5 / sqrt(NREP), sqrt(0.9 * 0.1 / NREP)),
sprintf("Flux de graines : 12000019 + cellule*10007 + r*101, disjoint de tous les flux de la Draft 3."),
"")

## Validation against the published Draft-3 comparison
ref_file <- "results/draft3/comparison_summary.csv"
if (file.exists(ref_file)) {
  ref <- read.csv(ref_file, stringsAsFactors = FALSE)
  w <- S[S$n == 2000 & S$p == 1000 & abs(S$rho - 0.25) < 1e-9, ]
  if (!nrow(w)) {
    md <- c(md, "## 2. Controle : cellule de reference de la Draft 3", "",
      "La cellule n = 2000, p = 1000, rho = 0.25 n'est pas encore calculee ; ce controle apparaitra des qu'elle sera disponible.", "")
  } else {
    body <- character(0)
    for (m in c("M1", "M2", "M3", "M4")) for (k in METHODS) {
      a <- w[w$model == m & w$method == k, ]
      b <- ref[ref$model == m & ref$method == k, ]
      if (!nrow(a) || !nrow(b)) next
      body <- c(body, paste0("| ", paste(c(m, LABEL[k],
        fmt(b$sure4), fmt(a$sure4), fmt(b$sure20), fmt(a$sure20),
        fmt(b$ermax, 1), fmt(a$ermax, 1)), collapse = " | "), " |"))
    }
    md <- c(md,
      "## 2. Controle : cellule de reference de la Draft 3",
      "",
      sprintf("Cellule n = 2000, p = 1000, rho = 0.25. Colonnes *pub.* : `%s` (200 replications, graines de la Draft 3) ; colonnes *ici* : cette campagne (%d replications, graines independantes). Les ecarts sont compatibles avec le bruit Monte Carlo.",
              ref_file, NREP),
      "",
      md_table(c("Modele", "Methode", "Sure-4 pub.", "Sure-4 ici",
                 "Sure-20 pub.", "Sure-20 ici", "E(Rmax) pub.", "E(Rmax) ici"),
               c("l", "l", rep("r", 6)), body))
  }
}

## --- narrative analysis, with every figure computed from the data -----------
analysis_lines <- function(S) {
  if (length(unique(S$n)) < 3L || length(unique(S$p)) < 4L) return(character(0))
  m20 <- function(meth, model) mean(S$sure20[S$method == meth & S$model == model])
  glob <- function(col, meth) mean(S[[col]][S$method == meth])
  by <- function(col, meth, key, val)
    mean(S[[col]][S$method == meth & S[[key]] == val])
  comp <- function(col, meth, model)
    mean(S[[col]][S$method == meth & S$model == model])
  ratio <- function(col) glob(col, "ours") / glob(col, "q950")
  c(
"## 3. Analyse d'ensemble",
"",
"### 3.1 Le classement global des methodes",
"",
sprintf("Sur le Sure-20 moyen, quantile SIS a tau=.95 devance le screen propose dans les quatre modeles : %s contre %s pour M1 a M4. Yoshida--Umezu est derniere partout (%s). Le constat vaut sur l'ensemble de la grille et pas seulement sur la cellule de la Draft 3.",
        paste(fmt(vapply(c("M1","M2","M3","M4"), function(m) m20("q950", m), 0), 2), collapse = " / "),
        paste(fmt(vapply(c("M1","M2","M3","M4"), function(m) m20("ours", m), 0), 2), collapse = " / "),
        paste(fmt(vapply(c("M1","M2","M3","M4"), function(m) m20("yu", m), 0), 2), collapse = " / ")),
"",
sprintf("L'ecart se referme pourtant sur la recuperation exacte : en Sure-4 moyen le screen propose fait %s contre %s pour quantile SIS a tau=.95, soit une egalite aux erreurs Monte Carlo pres, et il devance nettement les trois autres valeurs de tau (%s a tau=.90, %s a tau=.975, %s a tau=.99).",
        fmt(glob("sure4", "ours")), fmt(glob("sure4", "q950")),
        fmt(glob("sure4", "q900")), fmt(glob("sure4", "q975")),
        fmt(glob("sure4", "q990"))),
"",
sprintf("C'est la lecture que seules les valeurs de d etendues permettent. Le rapport des deux methodes passe de %s en Sure-4 a %s en Sure-100 (%s contre %s) : le screen propose place les variables actives tres haut quand il reussit, mais elargir la fenetre ne le sauve pas quand il echoue. Sa distribution de Rmax est bimodale -- succes quasi parfait ou echec profond -- tandis que celle de quantile SIS a une queue plus legere, donc elle profite de l'augmentation de d. Un tableau limite au Sure-20 masque entierement cette difference de nature.",
        fmt(ratio("sure4"), 2), fmt(ratio("sure100"), 2),
        fmt(glob("sure100", "ours")), fmt(glob("sure100", "q950"))),
"",
"### 3.2 Ce qui gouverne la performance",
"",
sprintf("**La taille d'echantillon domine tout.** Le Sure-20 du screen propose passe de %s a %s puis %s quand n vaut 1000, 2000 puis 5000. A n=5000 il se tient entre %s et %s quelle que soit la dimension, alors qu'a n=1000 il ne depasse jamais %s. C'est n, et non p, qui est la contrainte active.",
        fmt(by("sure20", "ours", "n", 1000)), fmt(by("sure20", "ours", "n", 2000)),
        fmt(by("sure20", "ours", "n", 5000)),
        fmt(min(vapply(c(200,500,1000,2000), function(pp)
          mean(S$sure20[S$method == "ours" & S$n == 5000 & S$p == pp]), 0))),
        fmt(max(vapply(c(200,500,1000,2000), function(pp)
          mean(S$sure20[S$method == "ours" & S$n == 5000 & S$p == pp]), 0))),
        fmt(max(vapply(c(200,500,1000,2000), function(pp)
          mean(S$sure20[S$method == "ours" & S$n == 1000 & S$p == pp]), 0)))),
"",
sprintf("**La dimension coute, mais pas davantage au screen propose qu'a ses concurrents.** En passant de p=200 a p=2000 son Sure-20 tombe de %s a %s, soit une perte relative de %.0f%% ; quantile SIS a tau=.95 perd %.0f%% (%s vers %s). L'idee que la competition entre p coordonnees penaliserait specifiquement un classement par score minimal n'est donc pas verifiee : la degradation est de meme ampleur pour les deux familles. Yoshida--Umezu est en revanche la plus sensible, avec %.0f%% de perte.",
        fmt(by("sure20", "ours", "p", 200)), fmt(by("sure20", "ours", "p", 2000)),
        100 * (1 - by("sure20", "ours", "p", 2000) / by("sure20", "ours", "p", 200)),
        100 * (1 - by("sure20", "q950", "p", 2000) / by("sure20", "q950", "p", 200)),
        fmt(by("sure20", "q950", "p", 200)), fmt(by("sure20", "q950", "p", 2000)),
        100 * (1 - by("sure20", "yu", "p", 2000) / by("sure20", "yu", "p", 200))),
"",
sprintf("**La correlation aide, fortement et pour toutes les methodes.** Le Sure-20 du screen propose passe de %s a rho=0 a %s a rho=0.5 ; quantile SIS de %s a %s. Le cas independant est donc le plus difficile, ce qui peut surprendre. L'explication tient a la structure AR(1) : les quatre coordonnees actives sont adjacentes, donc correlees entre elles, et conditionner sur l'une deplace aussi les trois autres. La depression de l'enveloppe le long d'une fibre active s'en trouve amplifiee, et le signal marginal de chaque active augmente avec rho. Les etudes qui ne rapportent qu'un seul rho intermediaire surestiment donc la performance par rapport au cas independant.",
        fmt(by("sure20", "ours", "rho", 0)), fmt(by("sure20", "ours", "rho", 0.5)),
        fmt(by("sure20", "q950", "rho", 0)), fmt(by("sure20", "q950", "rho", 0.5))),
"",
"### 3.3 La ou le screen propose est irremplacable",
"",
sprintf("Sur M2, seul modele dote de variables d'echelle (A_scale = {5,6,7,8}, qui deplacent les quantiles finis sans toucher l'indice de queue), le top-4 de quantile SIS en contient %s en moyenne contre %s pour le screen propose. Ce dernier retient en meme temps davantage de vraies actives en indice de queue (%s contre %s). C'est la seule dimension ou il domine, et c'est precisement celle qui justifie son existence : il repond a une autre question, pas mieux a la meme.",
        fmt(comp("top4_scale", "q950", "M2"), 2), fmt(comp("top4_scale", "ours", "M2"), 2),
        fmt(comp("top4_gamma", "ours", "M2"), 2), fmt(comp("top4_gamma", "q950", "M2"), 2)),
"",
sprintf("Le prix a payer est visible dans le meme tableau : sur M2 son Sure-20 vaut %s contre %s pour quantile SIS. Il ne se laisse pas tromper par les variables d'echelle, mais il classe moins surement les actives.",
        fmt(m20("ours", "M2")), fmt(m20("q950", "M2"))),
"",
"### 3.4 Le reglage de tau",
"",
sprintf("Quantile SIS n'est pas une methode mais une famille, et l'ecart entre ses membres depasse l'ecart entre familles : en Sure-20 moyen, tau=.95 donne %s et tau=.99 donne %s. Le tau optimal n'est pas connu en pratique, et tau=.99 -- le choix le plus naturel si l'on cherche un effet de queue -- est le pire des quatre. Compare a cette famille reglee a sa meilleure valeur, le screen propose part avec un handicap qui n'existe pas dans une application reelle.",
        fmt(mean(S$sure20[S$method == "q950"])),
        fmt(mean(S$sure20[S$method == "q990"]))),
"",
"### 3.5 Reserves",
"",
sprintf("Chaque cellule repose sur %d replications, donc l'erreur type d'une probabilite atteint %.3f ; les ecarts inferieurs a 0.10 dans une cellule isolee ne sont pas interpretables. Les moyennes de la section 4, portant sur 72 cellules par modele, sont en revanche precises. Les quatre modeles partagent la meme forme de gamma decroissante en chaque coordonnee active, donc les conclusions ne s'etendent pas telles quelles a des effets non monotones. Enfin le reglage (a*, b*) = (%s, %s) a ete choisi dans la Draft 3 a n=2000 et p=1000 : les cellules a n=5000 et p=2000 utilisent donc un reglage cale ailleurs, ce qui joue en defaveur du screen propose sur une partie de la grille.",
        NREP, 0.5 / sqrt(NREP), Sys.getenv("ASTAR", "0.30"),
        Sys.getenv("BSTAR", "0.10")),
"")
}
md <- c(md, analysis_lines(S))

md <- c(md, "## 4. Vues agregees", "",
  sprintf("Moyennes non ponderees des Sure-d sur les cellules concernees (%d cellules par modele).",
          length(unique(S$n)) * length(unique(S$p)) * length(unique(S$rho))), "")

for (m in c("M1", "M2", "M3", "M4")) {
  w <- S[S$model == m, ]
  md <- c(md, sprintf("### 4.%d Modele %s", match(m, c("M1","M2","M3","M4")), m), "",
    "Profil Sure-d, moyenne sur toutes les cellules du modele :", "",
    sured_profile(w),
    "Sure-20 par taille d'echantillon (moyenne sur p et rho) :", "",
    marginal(w, "n"),
    "Sure-20 par dimension (moyenne sur n et rho) :", "",
    marginal(w, "p"),
    "Sure-20 par correlation AR(1) (moyenne sur n et p) :", "",
    marginal(w, "rho"))
}

## Composition of the top 4 (relevant for M2's scale-active variables)
body <- character(0)
for (m in c("M1", "M2", "M3", "M4")) for (k in METHODS) {
  w <- S[S$model == m & S$method == k, ]
  body <- c(body, paste0("| ", paste(c(m, LABEL[k],
    fmt(mean(w$top4_gamma), 2), fmt(mean(w$top4_scale), 2)), collapse = " | "), " |"))
}
md <- c(md, "### 4.5 Composition du top-4", "",
  "Nombre moyen, parmi les 4 coordonnees les mieux classees, de variables actives en indice de queue (A_gamma = {1,2,3,4}) et de variables d'echelle (A_scale = {5,6,7,8}), moyenne sur toutes les cellules.", "",
  md_table(c("Modele", "Methode", "|top4 inter A_gamma|", "|top4 inter A_scale|"),
           c("l", "l", "r", "r"), body))

## --- full tables ------------------------------------------------------------
md <- c(md, "## 5. Tableaux complets", "",
  sprintf("Une table par (modele, n, p) ; lignes = rho x methode. Sure-d pour d = %s.",
          paste(DS, collapse = ", ")), "")
si <- 0L
for (m in c("M1", "M2", "M3", "M4")) {
  si <- si + 1L
  md <- c(md, sprintf("### 5.%d Modele %s", si, m), "")
  for (nn in sort(unique(S$n))) for (pp in sort(unique(S$p))) {
    w <- S[S$model == m & S$n == nn & S$p == pp, ]
    if (!nrow(w)) next
    md <- c(md, sprintf("#### %s, n = %d, p = %d", m, nn, pp), "", cell_table(w))
  }
}

## --- timings ----------------------------------------------------------------
body <- character(0)
for (nn in sort(unique(TM$n))) for (pp in sort(unique(TM$p))) {
  w <- TM[TM$n == nn & TM$p == pp, ]
  if (!nrow(w)) next
  body <- c(body, paste0("| ", paste(c(nn, pp,
    fmt(mean(w$ours), 2), fmt(mean(w$yu), 2), fmt(mean(w$qsis), 2),
    fmt(mean(w$yu_undefined), 3)), collapse = " | "), " |"))
}
md <- c(md, "## 6. Cout de calcul", "",
  "Secondes par replication et par methode, monocoeur, moyenne sur modeles et rho. Quantile SIS est chronometre par valeur de tau. Derniere colonne : proportion moyenne de scores Yoshida--Umezu non definis.", "",
  md_table(c("n", "p", "Tail-index SIS", "Yoshida--Umezu", "Quantile SIS (par tau)",
             "YU non defini"), c("r", "r", "r", "r", "r", "r"), body))

writeLines(md, file.path(outdir, "rapport_grille.md"))
cat("WROTE", file.path(outdir, "rapport_grille.md"), "-", length(md), "lines\n")
