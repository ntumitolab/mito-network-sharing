#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are two arguments: if not, return an error
if (length(args)!=2) {
  stop("Exactly two arguments must be supplied (input file).\n", call.=FALSE)
}

graphs.file = args[1]
results.file = args[2]
output.file = paste(args[1], ".prediction.png", sep="")

library(igraph)
library(ggplot2)

# graphs.file = "orig-1.xml-amlist.csv-0-0-graphs.txt"
# results.file = "orig-1.xml-amlist.csv-0-0-results-frame.txt"
graphs.df = read.table(graphs.file)
results.df = read.csv(results.file, header=T)

points.df = data.frame(L=NULL, expt=NULL, obs=NULL, est=NULL)
expts = unique(graphs.df$V1)
for(L in 2:10) {
  for(i in 1:length(expts)) {
    amlist = graphs.df[graphs.df$V1 == expts[i] & graphs.df$V2 == 0, 3:4]
    g = graph_from_edgelist(as.matrix(amlist), directed=F)
    bingos = results.df$num.bingos[results.df$L == L & results.df$m == 0 & results.df$n.edges == 1000 & results.df$expt == expts[i] & results.df$it == 0]
 
    # g is a igraph object containing our graph
    deg=degree(g)
   
    # aaa stores the frequencies of nodes with degree
    # higher than L-1
    aaa=table(deg[which(deg[]>=L-1)])
    pr=0
    if (length(aaa)==0) {
      estimatedBingos=0
      next
    }
    for (ii in 1:length(aaa)) {
      n=as.numeric(rownames(aaa)[ii])+1
      comb=L
      spr=1
      pl=1
      for (j in 1:L) {
        pl=pl*(-1)
        spr=spr+pl*dim(combn(comb,j))[2]*((comb-j)/comb)^n
      }
      pr=pr+aaa[[ii]]*spr
    }
    estimatedBingos=pr

    points.df = rbind(points.df, data.frame(L=L, expt=i, obs = bingos, est = estimatedBingos))
  }
}

png(output.file, width=800, height=800)

ggplot(data = points.df, aes(x = est, y = obs, color = factor(L))) +
  geom_point(size=3) +
  geom_abline(slope = 1, intercept = 0) +
  coord_trans(x="sqrt", y="sqrt") +
  xlab("Predicted bingos") + ylab("Observed bingos") + labs(colour = "L") +
  theme(axis.title.x = element_text(size=24), axis.title.y = element_text(size=24),
        axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
        legend.title = element_text(size=24), legend.text = element_text(size=18))  


dev.off()
