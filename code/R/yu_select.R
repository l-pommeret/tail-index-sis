## Selection of the Yoshida--Umezu tuning pair (k, h) from a tuning table.
##
## The competitor is given a PER-MODEL oracle-style tuning: for each model
## we keep the grid cell that maximises its own Sure-20, ties broken by the
## smaller mean worst active rank.  This is deliberately more generous than
## the treatment of the proposed screen, which is reported at a single
## (a*, b*) common to all four models.  The asymmetry is stated in the
## manuscript; it makes the comparison conservative in our disfavour.
##
## `tab` must have columns k, h, sure20, ermax (one row per grid cell).
yu_select_best <- function(tab) {
  stopifnot(all(c("k", "h", "sure20", "ermax") %in% names(tab)))
  o <- order(-tab$sure20, tab$ermax, tab$h, tab$k)
  tab[o[1L], , drop = FALSE]
}

## Same, applied per model; returns one row per model.
yu_select_by_model <- function(tab) {
  stopifnot("model" %in% names(tab))
  do.call(rbind, lapply(split(tab, tab$model), yu_select_best))
}
