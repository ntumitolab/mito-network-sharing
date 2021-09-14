#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("Exactly one argument must be supplied (input file).n", call.=FALSE)
}

input.fname = args[1];
output.fname = paste(args[1], ".graphs.png", sep="")

# graph analysis
library(igraph)
library(brainGraph)


all.titles = c("Bio", "Low diff", "Mid diff", "High diff", "Low diff+cyt 1", "Mid diff+cyt 1", "High diff+cyt 1", "Low diff+cyt 2", "Mid diff+cyt 2", "High diff+cyt 2", "Diff+cyt+inactive 1", "Diff+cyt+inactive 2", "Diff+cyt+inactive 3", "ER", "SF 1", "SF 2", "SF 3", "WS 1", "WS 2", "WS 3", "Geometric", "Star")
for(cs in (3+seq(from=0,to=7)*5)) {
  for(ct in 1:2) {
    all.titles = c(all.titles, paste("clique ", cs, "-", ct, sep=""))
    }
}

to.plot = c(0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 25)+1

titles = all.titles[to.plot]

n = length(to.plot)

png(output.fname, width=1600, height=1600)
par(mfrow=c(5,4))
results = data.frame(expt=NULL, iteration=NULL, n=NULL, e=NULL, efficiency=NULL, groups=NULL)
all.graphs = read.table(input.fname)
for(i in 1:n)
{
  graphs.by.expt = all.graphs[all.graphs[,1] == to.plot[i]-1 & all.graphs[,2] == 0, 3:4]
  n.node = graphs.by.expt[1,3]
  g = graph_from_edgelist(as.matrix(graphs.by.expt[-1,]), directed=F)
  g = add_vertices(g, n.node-gorder(g))
  drange = max(degree(g))-min(degree(g))
  groups = walktrap.community(g)
  ngroups = max(groups$membership)
  results = data.frame(expt=i, n=gorder(g), e=gsize(g), efficiency=efficiency(g, type="global"), groups=ngroups, cc = count_components(g))
#  plot(g, vertex.size=0, vertex.label=NA, main=paste(c(titles[i], "\n", "n = ", results$n, " e = ", results$e, " cc = ", results$cc, "\nnu = ", format(results$efficiency, digits=3), " ng = ", results$groups), collapse=""))
  plot(g, vertex.size=0, vertex.label=NA)
  #title(paste(c(titles[i], "\n", "n = ", results$n, " e = ", results$e, " cc = ", results$cc, "\nν = ", format(results$efficiency, digits=3)), collapse=""), cex.main=2)
  title(paste(c(titles[i], "\nν = ", format(results$efficiency, digits=3), " deg range = ", drange), collapse=""), cex.main=2)
}
dev.off()


output.fname = paste(args[1], ".degrees.png", sep="")
png(output.fname, width=1600, height=1600)
par(mfrow=c(5,4))
for(i in 1:n)
{
  graphs.by.expt = all.graphs[all.graphs[,1] == to.plot[i]-1 & all.graphs[,2] == 0, 3:4]
  n.node = graphs.by.expt[1,3]
  g = graph_from_edgelist(as.matrix(graphs.by.expt[-1,]), directed=F)
  g = add_vertices(g, n.node-gorder(g))
  groups = walktrap.community(g)
  ngroups = max(groups$membership)
  results = data.frame(expt=i, n=gorder(g), e=gsize(g), efficiency=efficiency(g, type="global"), groups=ngroups, cc = count_components(g))
  hist(degree(g), breaks=(seq(1:(max(degree(g))+2.5))-1.5), main=titles[i], cex.main=2, cex.lab=2, cex.axis=2, xlim=c(0,max(degree(g))+1))
  # hist(degree(g), breaks=(seq(1:(max(degree(g))+2.5))-1.5), main=paste(c(titles[i], "\n", "n = ", results$n, " e = ", results$e, " cc = ", results$cc, "\nnu = ", format(results$efficiency, digits=3), " ng = ", results$groups), collapse=""), xlim=c(0,max(degree(g))+1))
#  hist(degree(g), breaks=(seq(1:100)-1.5), main=paste(c(titles[i], "\n", "n = ", results$n, " e = ", results$e, " cc = ", results$cc, "\nnu = ", format(results$efficiency, digits=3), " ng = ", results$groups), collapse=""), xlim=c(0,40))
}
dev.off()