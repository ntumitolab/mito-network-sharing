#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("Exactly two arguments must be supplied (input file and m).n", call.=FALSE)
}

input.fname = args[1];
m = as.numeric(args[2]);
output.fname = paste(args[1], ".results-m-", m, ".png", sep="")

library(ggplot2)

all.titles = c("Bio", "Low diff", "Mid diff", "High diff", "Low diff+cyt 1", "Mid diff+cyt 1", "High diff+cyt 1", "Low diff+cyt 2", "Mid diff+cyt 2", "High diff+cyt 2", "Diff+cyt+slow 1", "Diff+cyt+slow 2", "Diff+cyt+slow 3", "ER", "SF 1", "SF 2", "SF 3", "WS 1", "WS 2", "WS 3", "Geometric", "Star")
for(cs in (3+seq(from=0,to=7)*5)) {
  for(ct in 1:2) {
    all.titles = c(all.titles, paste("clique ", cs, "-", ct, sep=""))
    }
}

results <- read.csv(input.fname, header=T)
results = results[results$m==m,]
bio.set = c(0,results$bingo.5.mean[results$expt==0])

results$relative.bingo.score = results$bingo.5.mean / bio.set[results$L]

png(output.fname, width=800,height=800)
#borders = c(2-0.5, 10+0.5, 12-0.5, 14+0.5, 18-0.5, 18+0.5)-1
borders = c(   which(all.titles=="ER")-0.5, which(all.titles=="ER")+0.5,
               which(all.titles=="WS 1")-0.5, which(all.titles=="WS 3")+0.5,
	       which(all.titles=="Star")-0.5, which(all.titles=="Star")+0.5    )-1

ggplot() +
  geom_rect(data = data.frame(xmin = borders[1], xmax = borders[2], ymin = 0, ymax = Inf), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="white", alpha = 0.9) +
  geom_rect(data = data.frame(xmin = borders[3], xmax = borders[4], ymin = 0, ymax = Inf), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="white", alpha = 0.9) +
  geom_rect(data = data.frame(xmin = borders[5], xmax = borders[6], ymin = 0, ymax = Inf), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="white", alpha = 0.9) +
  geom_rect(data = data.frame(xmin = -0.5, xmax = 0.5, ymin = 0, ymax = Inf), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha = 0.15) +
  geom_point(data = results, aes(x=expt, y=relative.bingo.score, group=factor(L), colour=factor(L))) +
  geom_line(data = results, aes(x=expt, y=relative.bingo.score, group=factor(L), colour=factor(L))) +
  geom_hline(yintercept=1) +
  scale_y_continuous(trans='log', breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)) +
  scale_x_continuous(breaks=seq(from=0,to=length(all.titles)-1), labels=all.titles) +
  xlab("Network") + ylab("Bingo score relative to bio network") + labs(colour = "L") +
  theme(axis.text.x = element_text(size=14, angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24), axis.text.y = element_text(size=18),
	legend.title = element_text(size=24), legend.text = element_text(size=18))  

mean.relative = rep(0, 49)
for(i in seq(1:50)) {
  mean.relative[i] = mean(log(results$relative.bingo.score[results$expt==(i-1)]))
}

dev.off()