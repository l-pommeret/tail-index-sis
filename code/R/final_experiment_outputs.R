args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("usage: final_experiment_outputs.R RESULTS_DIR OUTPUT_DIR")
rdir<-args[1]; odir<-args[2]; dir.create(odir,recursive=TRUE,showWarnings=FALSE)
files<-c("main_dimension.rds","dependence.rds","models.rds","scaling.rds","sensitivity.rds")
if(file.exists(file.path(rdir,"extension_p10000.rds"))) files<-c(files,"extension_p10000.rds")

# Replicate-level metrics, including all active ranks and score summaries.
extract_rows<-function(jobs) lapply(jobs,function(z) {
 ar<-rep(NA_integer_,4); if(length(z$active_ranks)) ar[seq_along(z$active_ranks)]<-z$active_ranks
 inactive<-setdiff(seq_len(z$p),z$active)
 data.frame(model=z$model,signal=z$signal,n=z$n,p=z$p,rho=z$rho,a=z$a,b=z$b,
 replicate=z$replicate,seed=z$seed,rank1=ar[1],rank2=ar[2],rank3=ar[3],rank4=ar[4],
 max_active_rank=if(length(z$active_ranks))max(z$active_ranks) else NA,
 active_score=if(length(z$active))mean(z$scores[z$active]) else NA,
 inactive_score=mean(z$scores[inactive]),under_rate=mean(z$under_rate),elapsed=z$elapsed)
})
raw_parts<-lapply(files,function(f) { x<-readRDS(file.path(rdir,f))$jobs; o<-do.call(rbind,extract_rows(x)); rm(x);gc();o })
raw<-do.call(rbind,raw_parts);rm(raw_parts);gc();write.csv(raw,file.path(odir,"replicate_metrics.csv"),row.names=FALSE)

# Special correlated-neighbour diagnostics for model A dependence experiment.
dep<-readRDS(file.path(rdir,"dependence.rds"))$jobs
nr<-do.call(rbind,lapply(dep,function(z) {
 sel<-z$ordering[seq_len(20)]; groups<-list(active=1:4,neighbour=5:8,remote=9:z$p)
 do.call(rbind,lapply(names(groups),function(g)data.frame(rho=z$rho,replicate=z$replicate,
 group=g,mean_score=mean(z$scores[groups[[g]]]),selected_top20=sum(groups[[g]]%in%sel))))
}))
key<-interaction(nr$rho,nr$group,drop=TRUE)
ns<-do.call(rbind,lapply(split(nr,key),function(w)data.frame(rho=w$rho[1],group=w$group[1],
 reps=length(unique(w$replicate)),mean_score=mean(w$mean_score),mean_selected_top20=mean(w$selected_top20))))
write.csv(ns,file.path(odir,"correlated_neighbours.csv"),row.names=FALSE)

# Null coordinate frequencies with Monte-Carlo standardized residuals.
null0<-readRDS(file.path(rdir,"models.rds"))$jobs
null<-null0[vapply(null0,function(z)z$model=="N",logical(1))]
if(length(null)) { p<-null[[1]]$p; dd<-20L; fr<-tabulate(unlist(lapply(null,function(z)z$ordering[1:dd])),nbins=p)
 ex<-length(null)*dd/p; no<-data.frame(coordinate=1:p,selected=fr,expected=ex,z=(fr-ex)/sqrt(ex*(1-dd/p)))
 write.csv(no,file.path(odir,"null_frequency_diagnostics.csv"),row.names=FALSE) }

# Publication-ready vector figures from audited raw records.
pdf(file.path(odir,"sure_screening_vs_p.pdf"),width=6.4,height=4.5)
md<-subset(raw,model=="A"&rho==0&a==.2&b==.2)
ag<-aggregate(md$max_active_rank<=20,list(n=md$n,p=md$p),mean); names(ag)[3]<-"sure"
plot(range(ag$p),c(0,1),type="n",log="x",xlab="Ambient dimension p",ylab="Sure-screening probability (d=20)")
for(nn in sort(unique(ag$n))) {w<-ag[ag$n==nn,];lines(w$p,w$sure,type="b",pch=match(nn,sort(unique(ag$n))))}
legend("topright",legend=paste0("n=",sort(unique(ag$n))),pch=seq_along(unique(ag$n)),lty=1,bty="n");dev.off()
pdf(file.path(odir,"sure_screening_vs_rho.pdf"),width=6.4,height=4.5)
dd<-subset(raw,model=="A"&n==1000&p==1000&a==.2&b==.2); ag<-aggregate(dd$max_active_rank<=20,list(rho=dd$rho),mean)
plot(ag$rho,ag$x,type="b",ylim=c(0,1),xlab=expression(rho),ylab="Sure-screening probability (d=20)",pch=19);dev.off()
pdf(file.path(odir,"sensitivity_heatmap.pdf"),width=7,height=5)
ss<-subset(raw,model=="A"&n==500&p==100&rho==0); ag<-aggregate(ss$max_active_rank<=20,list(a=ss$a,b=ss$b),mean)
image(sort(unique(ag$a)),sort(unique(ag$b)),matrix(ag$x,nrow=length(unique(ag$a))),xlab="tail exponent a",ylab="bandwidth exponent b",col=hcl.colors(20,"YlOrRd",rev=TRUE));dev.off()
cat("WROTE final tables and figures to",odir,"\n")
