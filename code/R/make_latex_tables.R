args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("usage: make_latex_tables.R TABLE_DIR MANUSCRIPT_TABLE_DIR")
td <- args[1]; od <- args[2]; dir.create(od, recursive=TRUE, showWarnings=FALSE)
fmt <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
se <- function(p,n) ifelse(is.na(p), NA_real_, sqrt(p*(1-p)/n))

md <- read.csv(file.path(td,"main_dimension.csv"))
md <- md[order(md$n,md$p),]
z <- c("\\begin{table}[t]","\\centering",
"\\caption{Model A over sample size and ambient dimension. Entries are probabilities with Monte Carlo standard errors in parentheses; each cell has 200 replications. Sure-20 means that all four active variables occur among the 20 smallest scores.}",
"\\label{tab:dimension}","\\begin{tabular}{rrrrr}","\\toprule",
"$n$ & $p$ & Sure-20 & Exact & Mean worst rank\\\\","\\midrule")
for(i in seq_len(nrow(md))) z <- c(z,sprintf("%d & %d & %s (%s) & %s (%s) & %.1f\\\\",
 md$n[i],md$p[i],fmt(md$sure20[i]),fmt(se(md$sure20[i],md$reps[i])),
 fmt(md$exact[i]),fmt(se(md$exact[i],md$reps[i])),md$mean_max_rank[i]))
z <- c(z,"\\bottomrule","\\end{tabular}","\\end{table}")
writeLines(z,file.path(od,"main_dimension.tex"))

dd <- read.csv(file.path(td,"dependence.csv")); dd<-dd[order(dd$rho),]
z <- c("\\begin{table}[t]","\\centering",
"\\caption{Model A with AR(1) covariates, $n=p=1000$. Monte Carlo standard errors are in parentheses; each cell has 200 replications.}",
"\\label{tab:dependence}","\\begin{tabular}{rrrr}","\\toprule",
"$\\rho$ & Sure-20 & Exact & Mean worst rank\\\\","\\midrule")
for(i in seq_len(nrow(dd))) z<-c(z,sprintf("%.2f & %s (%s) & %s (%s) & %.1f\\\\",dd$rho[i],
 fmt(dd$sure20[i]),fmt(se(dd$sure20[i],dd$reps[i])),fmt(dd$exact[i]),fmt(se(dd$exact[i],dd$reps[i])),dd$mean_max_rank[i]))
z<-c(z,"\\bottomrule","\\end{tabular}","\\end{table}");writeLines(z,file.path(od,"dependence.tex"))

mm <- read.csv(file.path(td,"models.csv")); mm<-mm[mm$model %in% c("B","C","D","E","F"),]
z <- c("\\begin{table}[t]","\\centering",
"\\caption{Tail-index designs at $n=p=1000$ and $\\rho=0$. Monte Carlo standard errors are in parentheses; each cell has 200 replications.}",
"\\label{tab:models}","\\begin{tabular}{lrrr}","\\toprule",
"Model & Sure-20 & Exact & Mean worst rank\\\\","\\midrule")
for(i in seq_len(nrow(mm))) { lab<-if(mm$model[i]=="B") sprintf("B ($c=%.2g$)",mm$signal[i]) else mm$model[i]; z<-c(z,sprintf("%s & %s (%s) & %s (%s) & %.1f\\\\",lab,
 fmt(mm$sure20[i]),fmt(se(mm$sure20[i],mm$reps[i])),fmt(mm$exact[i]),fmt(se(mm$exact[i],mm$reps[i])),mm$mean_max_rank[i]))
}
z<-c(z,"\\bottomrule","\\end{tabular}","\\end{table}");writeLines(z,file.path(od,"models.tex"))

## Screening probability as a function of the retained size d, for the designs
## whose exact-recovery probability vanishes at d = 20.
bd <- read.csv(file.path(td,"models_by_d.csv"))
bd <- bd[bd$n==1000 & bd$p==1000,]
sel <- list(c("B",1),c("B",0.5),c("C",1),c("D",1),c("E",1))
dgrid <- c(4,20,50,100,200,500)
lab <- c("A (=B, $c=1$)","B ($c=0.5$)","C","D","E")
z <- c("\\begin{table}[t]","\\centering",
"\\caption{Sure screening as a function of the retained size $d$, at $n=p=1000$ and $\\rho=0$ (200 replications). Upper block: probability that all active coordinates appear among the $d$ smallest scores. Lower block: mean number of active coordinates retained, to be compared with the chance value $ds_n/p$ shown in the last row of each design.}",
"\\label{tab:byd}","\\begin{tabular}{lrrrrrr}","\\toprule",
paste0("Design & ",paste(sprintf("$d=%d$",dgrid),collapse=" & "),"\\\\"),"\\midrule")
for(i in seq_along(sel)) { w<-bd[bd$model==sel[[i]][1] & bd$signal==as.numeric(sel[[i]][2]),]
  v<-sapply(dgrid,function(dd) w$sure[w$d==dd][1])
  z<-c(z,sprintf("%s & %s\\\\",lab[i],paste(fmt(v),collapse=" & "))) }
z <- c(z,"\\midrule","\\multicolumn{7}{l}{\\emph{Mean number of active coordinates retained}}\\\\")
for(i in seq_along(sel)) { w<-bd[bd$model==sel[[i]][1] & bd$signal==as.numeric(sel[[i]][2]),]
  v<-sapply(dgrid,function(dd) w$mean_true[w$d==dd][1])
  s<-length(unique(w$d)); sn <- if(sel[[i]][1]=="E") 2 else 4
  z<-c(z,sprintf("%s & %s\\\\",lab[i],paste(sprintf("%.2f",v),collapse=" & ")))
  if(i==length(sel)) z<-c(z,sprintf("\\emph{chance ($s_n=2$)} & %s\\\\",
     paste(sprintf("%.2f",2*dgrid/1000),collapse=" & "))) }
z <- c(z,sprintf("\\emph{chance ($s_n=4$)} & %s\\\\",paste(sprintf("%.2f",4*dgrid/1000),collapse=" & ")))
z <- c(z,"\\bottomrule","\\end{tabular}","\\end{table}")
writeLines(z,file.path(od,"by_d.tex"))

cc <- read.csv(file.path(td,"comparators.csv"))
z <- c("\\begin{table}[t]","\\centering",
"\\caption{Fixed-tuning common-data comparison at $n=1000,p=100$ (200 paired replications per design). Parentheses contain binomial Monte Carlo standard errors.}",
"\\label{tab:comparators}","\\begin{tabular}{lrlrr}","\\toprule",
"Model & $\\rho$ & Method & Sure-20 & Exact\\\\","\\midrule")
for(i in seq_len(nrow(cc))) z<-c(z,sprintf("%s & %.2f & %s & %s (%s) & %s (%s)\\\\",cc$model[i],cc$rho[i],cc$method[i],
 fmt(cc$sure20[i]),fmt(se(cc$sure20[i],200)),fmt(cc$exact[i]),fmt(se(cc$exact[i],200))))
z<-c(z,"\\bottomrule","\\end{tabular}","\\end{table}");writeLines(z,file.path(od,"comparators.tex"))
cat("WROTE LaTeX tables to",od,"\n")
