#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("Exactly two arguments must be supplied (input file and m).n", call.=FALSE)
}

input.fname = args[1];
m = as.numeric(args[2]);
output.fname = paste(args[1], ".trace-m-", m, ".png", sep="")

library(ggplot2)

results = read.csv(input.fname,header=T)
subset = results[results$m==m & results$expt==0,]
subset$L = factor(subset$L)
png(output.fname, width=800,height=800)
ggplot(data = subset, aes(x=prop.edges,y=prop.bingos,group=L)) +
  geom_point(aes(color=L)) +
  stat_smooth(aes(group=L, color=L, fill=L)) +
  xlab("Proportion of edges used") + ylab("Bingo score") + labs(colour = "L", fill="L") +
  theme(axis.title.x = element_text(size=24), axis.title.y = element_text(size=24),
        axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
        legend.title = element_text(size=24), legend.text = element_text(size=18))

dev.off()

