args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2L) stop("usage: run_comparators.R CONFIG.csv OUTPUT.rds [CORES]")
cores <- if(length(args)>=3L) as.integer(args[3]) else 1L
source("code/R/generate.R"); source("code/R/yoshida_umezu.R")
library(Rcpp); sourceCpp("code/src/local_hill.cpp", rebuild=FALSE)
cfg <- read.csv(args[1], stringsAsFactors=FALSE)
grid <- do.call(rbind,lapply(seq_len(nrow(cfg)),function(i)data.frame(i=i,rep=seq_len(cfg$reps[i]))))
one <- function(ix) {
  i<-grid$i[ix];rep<-grid$rep[ix]
  z <- cfg[i,]; seed <- 7000003L+i*10007L+rep*101L
  d <- simulate_dataset(z$n,z$p,z$rho,z$model,z$signal,seed,z$block_size)
  alpha <- z$n^(-z$a); h <- z$n^(-z$b)/2
  t0 <- proc.time()[3]
  ours <- vapply(seq_len(z$p),function(j) score_coordinate_cpp(d$x[,j],d$y,h,alpha,z$epsilon)[1],0.0)
  ours_time <- proc.time()[3]-t0
  k <- floor(z$yu_k_fraction*z$n)
  # Yoshida--Umezu is defined on the empirical marginal transform, matching
  # the r/(n+1) convention used elsewhere.  `d$u` is the population PIT used
  # only to generate the response and must never be passed to an estimator.
  yu_u <- rank_uniform_matrix(d$x)
  t0 <- proc.time()[3]; yu <- yu_score_matrix(yu_u,d$y,k,z$yu_h,
    seq(.02,.98,length.out=z$yu_grid)); yu_time <- proc.time()[3]-t0
  active <- if(z$model=="E") 1:2 else 1:4
  ord_o <- order(ours,seq_along(ours)); ord_y <- order(-yu$scores,seq_along(yu$scores),na.last=TRUE)
  list(config_row=i,replicate=rep,seed=seed,model=z$model,n=z$n,p=z$p,
    rho=z$rho,active=active,ours_scores=ours,yu_scores=yu$scores,
    ours_ordering=ord_o,yu_ordering=ord_y,ours_active_ranks=match(active,ord_o),
    yu_active_ranks=match(active,ord_y),yu_undefined=mean(!is.finite(yu$scores)),
    ours_elapsed=unname(ours_time),yu_elapsed=unname(yu_time),k=k,yu_h=z$yu_h,
    yu_empirical_pit=TRUE)
}
if(.Platform$OS.type=="unix"&&cores>1L) {
  jobs<-parallel::mclapply(seq_len(nrow(grid)),one,mc.cores=cores,mc.preschedule=FALSE)
} else {
  jobs<-lapply(seq_len(nrow(grid)),one)
}
dir.create(dirname(args[2]),recursive=TRUE,showWarnings=FALSE)
saveRDS(list(meta=list(config=cfg,r_version=R.version.string),jobs=jobs),args[2],compress="xz")
cat("WROTE",args[2],"JOBS",length(jobs),"\n")
