args <- commandArgs(trailingOnly=TRUE); stopifnot(length(args)==2L)
x <- readRDS(args[1])$jobs
key <- vapply(x,function(z) paste(z$model,z$n,z$p,z$rho,sep="|"),"")
out <- do.call(rbind,lapply(split(x,key),function(w) {
 z<-w[[1]]; one <- function(method) {
   ranks<-lapply(w,`[[`,paste0(method,"_active_ranks"))
   c(sure20=mean(vapply(ranks,max,0)<=20),exact=mean(vapply(ranks,max,0)<=length(z$active)),
     mean_max_rank=mean(vapply(ranks,max,0)),mean_elapsed=mean(vapply(w,`[[`,0.0,paste0(method,"_elapsed"))))
 }
 rbind(data.frame(model=z$model,n=z$n,p=z$p,rho=z$rho,method="Local-Hill SIS",t(one("ours"))),
       data.frame(model=z$model,n=z$n,p=z$p,rho=z$rho,method="Yoshida-Umezu",t(one("yu"))))
}))
write.csv(out,args[2],row.names=FALSE)
