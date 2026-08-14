## Draft-2 real-data application: Communities and Crime (Unnormalized), UCI.
## Response: violent crimes per 100k population (violentPerPop, column 146).
## Response-blind preprocessing:
##   rows: response observed and positive;
##   predictors: columns 6..129 of the raw file (the documented predictive
##     attributes), keeping numeric columns with < 5% missing values and at
##     least 100 distinct observed values (high cardinality, so empirical
##     ranks are meaningful).  No predictor is chosen using the response.
## Missing predictor values are handled pairwise: coordinate j is scored on
## the rows where X_j is observed.
## usage: Rscript code/R/real_crime.R ASTAR BSTAR OUTPUT.rds

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
ASTAR <- as.numeric(args[1L]); BSTAR <- as.numeric(args[2L])
out_path <- args[3L]

suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/qa_sis.R"); source("code/R/yoshida_umezu.R")

raw <- read.csv("data/crime/CommViolPredUnnormalizedData.txt",
                header = FALSE, na.strings = "?")
## Official ordered attribute list (UCI id 211): 5 identifiers, 124
## predictive attributes, 18 potential goals.  Hyphens sanitized.
attr_names <- c("communityname","state","countyCode","communityCode","fold",
  "pop","perHoush","pctBlack","pctWhite","pctAsian","pctHisp","pct12to21",
  "pct12to29","pct16to24","pct65up","persUrban","pctUrban","medIncome",
  "pctWwage","pctWfarm","pctWdiv","pctWsocsec","pctPubAsst","pctRetire",
  "medFamIncome","perCapInc","whitePerCap","blackPerCap","NAperCap",
  "asianPerCap","otherPerCap","hispPerCap","persPoverty","pctPoverty",
  "pctLowEdu","pctNotHSgrad","pctCollGrad","pctUnemploy","pctEmploy",
  "pctEmployMfg","pctEmployProfServ","pctOccupManu","pctOccupMgmt",
  "pctMaleDivorc","pctMaleNevMar","pctFemDivorc","pctAllDivorc","persPerFam",
  "pct2Par","pctKids2Par","pctKids4w2Par","pct12to17w2Par","pctWorkMom6",
  "pctWorkMom18","kidsBornNevrMarr","pctKidsBornNevrMarr","numForeignBorn",
  "pctFgnImmig3","pctFgnImmig5","pctFgnImmig8","pctFgnImmig10","pctImmig3",
  "pctImmig5","pctImmig8","pctImmig10","pctSpeakOnlyEng","pctNotSpeakEng",
  "pctLargHousFam","pctLargHous","persPerOccupHous","persPerOwnOccup",
  "persPerRenterOccup","pctPersOwnOccup","pctPopDenseHous","pctSmallHousUnits",
  "medNumBedrm","houseVacant","pctHousOccup","pctHousOwnerOccup",
  "pctVacantBoarded","pctVacant6up","medYrHousBuilt","pctHousWOphone",
  "pctHousWOplumb","ownHousLowQ","ownHousMed","ownHousUperQ","ownHousQrange",
  "rentLowQ","rentMed","rentUpperQ","rentQrange","medGrossRent",
  "medRentpctHousInc","medOwnCostpct","medOwnCostPctWO","persEmergShelt",
  "persHomeless","pctForeignBorn","pctBornStateResid","pctSameHouse5",
  "pctSameCounty5","pctSameState5","numPolice","policePerPop",
  "policeField","policeFieldPerPop","policeCalls","policCallPerPop",
  "policCallPerOffic","policePerPop2","racialMatch","pctPolicWhite",
  "pctPolicBlack","pctPolicHisp","pctPolicAsian","pctPolicMinority",
  "officDrugUnits","numDiffDrugsSeiz","policAveOT","landArea","popDensity",
  "pctUsePubTrans","policCarsAvail","policOperBudget","pctPolicPatrol",
  "gangUnit","pctOfficDrugUnit","policBudgetPerPop",
  "murders","murdPerPop","rapes","rapesPerPop","robberies","robbbPerPop",
  "assaults","assaultPerPop","burglaries","burglPerPop","larcenies",
  "larcPerPop","autoTheft","autoTheftPerPop","arsons","arsonsPerPop",
  "violentPerPop","nonViolPerPop")
stopifnot(length(attr_names) == ncol(raw))
names(raw) <- attr_names

y_all <- raw$violentPerPop
keep_row <- !is.na(y_all) & y_all > 0
pred_all <- raw[, 6:129]
miss_frac <- colSums(is.na(pred_all)) / nrow(pred_all)
n_distinct <- vapply(pred_all, function(v)
  length(unique(v[!is.na(v)])), integer(1))
sel <- miss_frac < 0.05 & n_distinct >= 100
X <- as.matrix(pred_all[keep_row, sel])
y <- y_all[keep_row]
n <- length(y); p <- ncol(X)
alpha <- n^(-ASTAR); h <- n^(-BSTAR) / 2

## Harmonized response-blind preprocessing: the (at most one) missing value
## per retained predictor is median-imputed, and every method (proposed
## screen and competitors) works from the same imputed matrix.
for (j in seq_len(p)) X[is.na(X[, j]), j] <- median(X[, j], na.rm = TRUE)

score_all <- function(a, b, rows = seq_len(n), tie_seed = NULL) {
  al <- length(rows)^(-a); hh <- length(rows)^(-b) / 2
  vapply(seq_len(p), function(j) {
    xj <- X[rows, j]
    if (!is.null(tie_seed)) {
      set.seed(tie_seed + j)
      xj <- xj + runif(length(xj), -1e-9, 1e-9) * (abs(xj) + 1)
    }
    score_coordinate_cpp(xj, y[rows], hh, al, 0.05)[["score"]]
  }, numeric(1))
}

## Primary fit and tuning-stability grid around (ASTAR, BSTAR).
primary <- score_all(ASTAR, BSTAR)
grid <- expand.grid(a = ASTAR + c(-.05, 0, .05), b = BSTAR + c(-.05, 0, .05))
grid_ranks <- lapply(seq_len(nrow(grid)), function(i) {
  s <- score_all(grid$a[i], grid$b[i]); rank(s, ties.method = "first")
})

## Random tie-breaking audit (predictors have a few exact ties).
tie_ranks <- lapply(1:20, function(s)
  rank(score_all(ASTAR, BSTAR, tie_seed = 900000 + s), ties.method = "first"))

## 200 seeded random half-samples.
half_ranks <- lapply(1:200, function(s) {
  set.seed(500000 + s)
  rows <- sort(sample(n, floor(n / 2)))
  rank(score_all(ASTAR, BSTAR, rows = rows), ties.method = "first")
})

uh <- apply(X, 2L, function(v) rank(v, ties.method = "average") / (n + 1))
yu <- yu_score_matrix(uh, y, k = floor(0.072 * n), h = 1)
q95 <- qa_sis_scores(uh, y, tau = 0.95)

## Hill diagnostic and Pareto QQ data for the response.
ys <- sort(y, decreasing = TRUE)
hill_k <- sort(unique(c(seq(20L, 600L, by = 5L), 204L)))
hill <- vapply(hill_k, function(k) mean(log(ys[1:k])) - log(ys[k + 1]),
               numeric(1))
kqq <- 600L
qq <- data.frame(x = log((n + 1) / (1:kqq)), y = log(ys[1:kqq]))

saveRDS(list(
  n = n, p = p, astar = ASTAR, bstar = BSTAR, alpha = alpha, h = h,
  nalpha = n * alpha, var_names = colnames(X),
  n_excluded_rows = sum(!keep_row),
  n_excluded_missing = sum(miss_frac >= 0.05),
  n_excluded_lowcard = sum(miss_frac < 0.05 & n_distinct < 100),
  primary_scores = primary, grid = grid, grid_ranks = grid_ranks,
  tie_ranks = tie_ranks, half_ranks = half_ranks,
  yu_scores = yu$scores, q95_scores = q95,
  hill_k = hill_k, hill = hill, qq = qq), out_path, compress = "xz")
cat("WROTE", out_path, ": n", n, "p", p, "nalpha", round(n * alpha), "\n")
