args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1L] else "results/yu_tuning_grid.rds"
z <- readRDS(path)
stopifnot(nrow(z$tuning) == 9L, nrow(z$raw) == 7200L,
          nrow(z$summary) == 36L,
          identical(sort(unique(z$raw$k_fraction)), c(.05, .072, .10)),
          identical(sort(unique(z$raw$h)), c(.5, 1, 1.5)),
          all(is.finite(as.matrix(z$raw[c("sure20", "exact", "mean_max_rank", "undefined")]))),
          all(z$raw$undefined == 0),
          all(table(interaction(z$raw$model, z$raw$rho, drop = TRUE),
                    z$raw$k_fraction, z$raw$h) == 200L))
cat("PASS Yoshida--Umezu tuning grid: 7,200 common-protocol fits, no undefined scores\n")
