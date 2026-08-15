## Ranking rules for the aggregated tail-index screen: audit and fix of the
## column-index tie-break.
##
## OLD rule (campaign5.R / campaign6.R / real_crime5.R):
##   per-setting ranks   rank(S[, k], ties.method = "first")   -> index ties
##   aggregate           amin = rowMin(ranks)
##   final order         order(amin, seq_len(p))               -> index ties
## The minimum of 9 rank permutations produces structural ties (up to 9
## coordinates can share amin = v), and with the true actives at columns
## 1--4 every active/inactive tie is resolved in favour of the active.
##
## NEW rule (permutation-invariant):
##   per-setting ranks   rank(S[, k], ties.method = "average")
##   lexicographic key   (min, median, mean, max) of the 9 ranks
##   final order         order(key..., tiebreak) where tiebreak is a
##                       seeded uniform draw, never the column index.
## The primary key is still the minimum rank: the aggregation philosophy
## (retain a coordinate as soon as one reasonable tuning supports it) is
## unchanged; the refinement only decides among coordinates with the SAME
## minimum, using their whole 9-setting rank profile, and residual exact
## ties (identical profiles) fall to a reproducible seeded draw.
##
## All screens in the study are marginal (the score of coordinate j is a
## function of (X_j, Y) alone), so score vectors are exactly equivariant
## under column permutation; the ranking rule is the only place where the
## column labels can leak into the results.

## p x K matrix of per-setting ranks, permutation-invariant.
setting_ranks_avg <- function(S)
  apply(S, 2L, rank, ties.method = "average", na.last = TRUE)

## Lexicographic aggregation key from a p x K rank matrix.
agg_key <- function(R)
  cbind(min  = apply(R, 1L, min),
        med  = apply(R, 1L, median),
        mean = rowMeans(R),
        max  = apply(R, 1L, max))

## Seeded uniform tie-break vector (isolated RNG state).
tiebreak_u <- function(p, seed) {
  old <- if (exists(".Random.seed", .GlobalEnv))
    get(".Random.seed", .GlobalEnv) else NULL
  set.seed(seed); u <- runif(p)
  if (is.null(old)) rm(".Random.seed", envir = .GlobalEnv) else
    assign(".Random.seed", old, .GlobalEnv)
  u
}

## --- OLD orderings (exact reproduction of campaign6.R lines 45-54) ------
order_selected_old <- function(S, idx_pub, p = nrow(S))
  order(S[, idx_pub], seq_len(p), na.last = TRUE)
order_agg_old <- function(S, p = nrow(S)) {
  R <- apply(S, 2L, rank, ties.method = "first")
  order(apply(R, 1L, min), seq_len(p))
}

## --- NEW orderings (permutation-invariant) ------------------------------
order_selected_new <- function(S, idx_pub, u)
  order(S[, idx_pub], u, na.last = TRUE)
order_agg_new <- function(S, u) {
  key <- agg_key(setting_ranks_avg(S))
  order(key[, 1L], key[, 2L], key[, 3L], key[, 4L], u)
}

## Bounds used by the tie diagnostic: same primary key as the OLD rule
## (first-ties min-rank), ties resolved for/against the active set.
order_agg_bound <- function(S, actives, favourable = TRUE, p = nrow(S)) {
  R <- apply(S, 2L, rank, ties.method = "first")
  amin <- apply(R, 1L, min)
  act <- as.integer(seq_len(p) %in% actives)
  if (favourable) order(amin, 1L - act, seq_len(p))
  else            order(amin,      act, seq_len(p))
}

## amin + pure random tie-break (no lexicographic refinement), to separate
## the effect of the refinement from the effect of dropping the index.
order_agg_minrandom <- function(S, u) {
  R <- apply(S, 2L, rank, ties.method = "first")
  order(apply(R, 1L, min), u)
}

## --- metrics ------------------------------------------------------------
metrics_from_order <- function(o, actives, scale_set)
  c(rmax = max(match(actives, o)),
    t4g  = sum(o[1:4]  %in% actives),
    t4s  = sum(o[1:4]  %in% scale_set),
    t24s = sum(o[1:24] %in% scale_set))

## --- tie diagnostics on one score matrix --------------------------------
## Returns counts of exact score ties, structural amin ties, and whether
## the tie group straddling positions 4 and 20 mixes active and inactive
## coordinates (the configurations in which the tie-break can move the
## reported metrics).
tie_stats <- function(S, actives) {
  p <- nrow(S)
  exact_score_ties <- sum(vapply(seq_len(ncol(S)), function(k)
    sum(duplicated(S[, k])), integer(1)))
  n_inf <- sum(!is.finite(S))
  R <- apply(S, 2L, rank, ties.method = "first")
  amin <- apply(R, 1L, min)
  tab <- table(amin)
  tied_coords <- sum(tab[tab > 1L])
  straddle <- function(d) {
    o <- order(amin, seq_len(p))
    v <- amin[o[d]]
    grp <- which(amin == v)
    inside  <- sum(amin < v)          # coords strictly better than the group
    cross   <- inside < d && inside + length(grp) > d
    mixed   <- any(grp %in% actives) && any(!(grp %in% actives))
    c(cross = cross, mixed_cross = cross && mixed, group_size = length(grp))
  }
  s4 <- straddle(4L); s20 <- straddle(20L)
  c(exact_score_ties = exact_score_ties, n_inf = n_inf,
    tied_coords_amin = tied_coords, n_distinct_amin = length(tab),
    cross4 = s4[["cross"]], mixed4 = s4[["mixed_cross"]],
    grp4 = s4[["group_size"]],
    cross20 = s20[["cross"]], mixed20 = s20[["mixed_cross"]],
    grp20 = s20[["group_size"]])
}
