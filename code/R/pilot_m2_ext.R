suppressMessages(library(Rcpp))
sourceCpp("code/src/local_hill.cpp", rebuild = FALSE)
source("code/R/generate3.R"); source("code/R/qa_sis.R")
n<-2000L;p<-1000L;rho<-0.25;h<-n^(-0.10)/2;alpha<-n^(-0.30)
ki<-6L; kap<-0.20
assign("KAPPA",kap,envir=globalenv())
out<-parallel::mclapply(1:80,function(r){
  seed<-7000003L+ki*100057L+r*503L
  d<-simulate_dataset3(n,p,rho,"M2",seed)
  sc<-vapply(seq_len(p),function(j) score_coordinate_cpp(d$z[,j],d$y,h,alpha,0.05),numeric(5))["score",]
  o_t<-order(sc,seq_len(p),na.last=TRUE)
  q95<-qa_sis_scores(apply(d$z,2L,function(x) rank(x,ties.method="average")/(n+1)),d$y,tau=0.95)
  o_q<-order(-q95,seq_len(p))
  list(kappa=kap,replicate=r,seed=seed,rt=match(1:4,o_t),rq=match(1:4,o_q),
       top4_t=o_t[1:4],top4_q=o_q[1:4])
},mc.cores=6,mc.preschedule=FALSE)
res<-readRDS("results/draft3/pilot_m2.rds"); res[["0.2"]]<-out
saveRDS(res,"results/draft3/pilot_m2.rds",compress="xz")
x<-out; rt<-sapply(x,function(z)max(z$rt)); rq<-sapply(x,function(z)max(z$rq))
cat(sprintf("kappa=0.20 Tail: S4=%.2f S20=%.2f Rmax=%.0f | Q95: S4=%.2f S20=%.2f Rmax=%.0f | top4 Tail g=%.2f s=%.2f | Q95 g=%.2f s=%.2f\n",
 mean(rt<=4),mean(rt<=20),mean(rt),mean(rq<=4),mean(rq<=20),mean(rq),
 mean(sapply(x,function(z)sum(z$top4_t%in%1:4))),mean(sapply(x,function(z)sum(z$top4_t%in%5:8))),
 mean(sapply(x,function(z)sum(z$top4_q%in%1:4))),mean(sapply(x,function(z)sum(z$top4_q%in%5:8)))))
