#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are two arguments: if not, return an error
if (length(args)!=2) {
  stop("Exactly two arguments must be supplied (input graphs and results files).\n", call.=FALSE)
}

#args = c("Output/orig-1.xml-amlist.csv-0-0-graphs.txt", "Output/orig-1.xml-amlist.csv-0-0-results-overall.txt")

graphs.file = args[1]
results.file = args[2]
output.file = paste(args[1], ".predicted.png", sep="")

library(igraph)
library(ggplot2)
library(gridExtra)

comparison= data.frame(observed=numeric(), 
                       estimated=numeric(),
                       threshold=numeric(),
                       L=numeric(),stringsAsFactors=F) 

# read the bingio results
bb=read.table(results.file,sep=",",header=T)
bb=bb[,c(1,2,3,30)]


for (L in 2:10) {
  
  results = data.frame(expt=NULL, iteration=NULL, n=NULL, e=NULL, efficiency=NULL, groups=NULL)
  all.graphs = read.table(graphs.file)
  graphs=list()
  numberOfBingos=vector()
  estimatedBingos=vector()
  
  for(i in 1:max(all.graphs[,1])+1)
  {
    graphs.by.expt = all.graphs[all.graphs[,1] == i-1 & all.graphs[,2] == 0, 3:4]
    n.node = graphs.by.expt[1,3]
    g = graph_from_edgelist(as.matrix(graphs.by.expt[-1,]), directed=F)
    g = add_vertices(g, n.node-gorder(g))
    graphs[[i]]=g
    numberOfBingos[i] = bb[bb[,1] == L & bb[,2] == 0 & bb[,3] == i-1, 4]
    
    # the estimated probability for the number of bingos
    deg=degree(g)
    aaa=table(deg[which(deg[]>=L-1)])
    pr=0
    if (length(aaa)==0) {
      estimatedBingos[i]=0
      next
    }
    for (ii in 1:length(aaa)) {
      n=as.numeric(rownames(aaa)[ii])+1
      comb=L
      spr=0
      pl=-1
      for (j in 0:L) {
        pl=pl*(-1)
        spr=spr+pl*dim(combn(comb,j))[2]*((comb-j)/comb)^n
      }
      pr=pr+aaa[[ii]]*spr
    }
    estimatedBingos[i]=pr
    
    
    # the approximate estimation using the expected value
    spr=0
    pl=-1
    for (j in 1:L) {
      pl=pl*(-1)
      spr=spr+pl*dim(combn(comb,j))[2]*1/(1-((comb-j)/comb))
    }
    spr2=0
    # spr=L*harmonic(L) #eventually the same with the above spr
    for (ii in 1:length(aaa)) {
      n=as.numeric(rownames(aaa)[ii])+1
      if (n>=spr) {
        spr2=spr2+aaa[[ii]]
      }
    }
    
    #store the estimated and the rough approximation (using the scalar approach)
    comparison[nrow(comparison) + 1,] = c(numberOfBingos[i],estimatedBingos[i],spr2,L)
  }
}


png(output.file, width=1000, height=500)

my.size.l = 18
my.size.m = 14
p1 = ggplot(data = comparison, aes(x = estimated, y = observed, color = factor(L))) +
  geom_point(size=3) +
  geom_abline(slope = 1, intercept = 0) + 
  xlab("Predicted bingos") + ylab("Observed bingos") + labs(colour = "L") +
  theme(axis.title.x = element_text(size=my.size.l), axis.title.y = element_text(size=my.size.l), axis.text.x = element_text(size=my.size.m), axis.text.y = element_text(size=my.size.m), legend.title = element_text(size=my.size.l), legend.text = element_text(size=my.size.m))
p2 = ggplot(data = comparison, aes(x = threshold, y = observed, color = factor(L))) +
  geom_point(size=3) +
  geom_abline(slope = 1, intercept = 0) + 
  xlab("Predicted bingos") + ylab("Observed bingos") + labs(colour = "L") +
  theme(axis.title.x = element_text(size=my.size.l), axis.title.y = element_text(size=my.size.l), axis.text.x = element_text(size=my.size.m), axis.text.y = element_text(size=my.size.m), legend.title = element_text(size=my.size.l), legend.text = element_text(size=my.size.m))
grid.arrange(p1, p2, nrow=1)

dev.off()
