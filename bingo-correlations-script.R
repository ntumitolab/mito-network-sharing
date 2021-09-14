#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("Exactly one argument must be supplied (input file).n", call.=FALSE)
}

#args = c("Output/orig-1.xml-amlist.csv-0-0-0-results-overall.txt")

input.fname = args[1];
output.fname = paste(args[1], ".correlations.png", sep="")

library(ggplot2)
library(GGally)

df = read.csv(input.fname, header=T)
subset = subset(df[(df$L==2 | df$L==5) & df$m==0,], select=c("L", "expt", "sd.degree.mean", "range.degree.mean", "efficiency.mean", "modularity.mean", "singleton.count.mean", "small.count.mean", "num.cc.mean", "mean.cc.size.mean", "bingo.1.mean", "bingo.3.mean", "bingo.5.mean"))

png(output.fname, width=800, height=800)
ggpairs(subset, mapping=ggplot2::aes(colour=factor(L)))
dev.off()
