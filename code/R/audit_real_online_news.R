args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1L] else "results/real_online_news/online_news_application.rds"
z <- readRDS(path)
stopifnot(is.list(z), nrow(z$scores) == 9L * 15L,
          nrow(z$split_scores) == 2L * 15L,
          nrow(z$tie_scores) == 20L * 15L,
          nrow(z$half_scores) == 200L * 15L,
          length(z$metadata$features) == 15L,
          z$metadata$analysis_rows == 38462L,
          all(is.finite(z$scores$score)), all(is.finite(z$split_scores$score)),
          all(is.finite(z$tie_scores$score)),
          all(is.finite(z$half_scores$score)),
          all(z$scores$under_rate == 0),
          all(vapply(split(z$scores$rank, interaction(z$scores$a, z$scores$b)),
                     function(v) identical(sort(v), seq_len(15L)), logical(1))),
          all(vapply(split(z$tie_scores$rank, z$tie_scores$seed),
                     function(v) identical(sort(v), seq_len(15L)), logical(1))),
          all(vapply(split(z$half_scores$rank, z$half_scores$seed),
                     function(v) identical(sort(v), seq_len(15L)), logical(1))))
cat("PASS real online-news audit: 38462 rows, 15 features, 9 tunings,",
    "2 temporal halves, 20 tie-break audits, 200 random half-samples\n")
