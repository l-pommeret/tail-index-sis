args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("usage: audit_null.R INPUT.rds OUTPUT.csv")
x <- readRDS(args[1L])$jobs
x <- x[vapply(x,function(z) z$model=="N",logical(1))]
stopifnot(length(x)>0)
p <- x[[1]]$p; d <- 20L
freq <- tabulate(unlist(lapply(x,function(z) z$ordering[seq_len(d)])),nbins=p)
expected <- length(x)*d/p
chisq <- sum((freq-expected)^2/expected)
pvalue <- pchisq(chisq,df=p-1,lower.tail=FALSE)
out <- data.frame(coordinate=seq_len(p),selected=freq,expected=expected)
write.csv(out,args[2L],row.names=FALSE)
cat("NULL reps",length(x),"p",p,"d",d,"range",range(freq),
    "chi_square",chisq,"df",p-1,"p_value",pvalue,"\n")
