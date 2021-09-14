#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("Exactly one argument must be supplied (input file).n", call.=FALSE)
}

input.fname = args[1];
output.fname = paste(args[1], ".graph.png", sep="")

# graph analysis
library(igraph)
library(brainGraph)


png(output.fname, width=800, height=800)
results = data.frame(expt=NULL, iteration=NULL, n=NULL, e=NULL, efficiency=NULL, groups=NULL)
all.graphs = read.table(input.fname)

graphs.by.expt = all.graphs[all.graphs[,1] == 0 & all.graphs[,2] == 0, 3:4]
n.node = graphs.by.expt[1,3]
g = graph_from_edgelist(as.matrix(graphs.by.expt[-1,]), directed=F)
g = add_vertices(g, n.node-gorder(g))
groups = walktrap.community(g)
ngroups = max(groups$membership)
results = data.frame(n=gorder(g), e=gsize(g), efficiency=efficiency(g, type="global"), groups=ngroups, cc = count_components(g))
plot(g, vertex.size=0, vertex.label=NA)
title(paste(c("n = ", results$n, " e = ", results$e, " cc = ", results$cc), collapse=""), cex.main=2)
dev.off()


