#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("Exactly two argument must be supplied (input file and m).n", call.=FALSE)
}

input.fname = args[1];
m = as.numeric(args[2]);
output.fname = paste(args[1], ".traces-m-", m, ".png", sep="")

library(ggplot2)

results = read.csv(input.fname,header=T)
subset = results[results$m==m & results$prop.bingos < 1,]
all.titles = c("X Bio", "Sim: Low diff", "Sim: Mid diff", "Sim: High diff", "Sim: Low diff+cyt 1", "Sim: Mid diff+cyt 1", "Sim: High diff+cyt 1", "Sim: Low diff+cyt 2", "Sim: Mid diff+cyt 2", "Sim: High diff+cyt 2", "Sim: Diff+cyt+slow 1", "Sim: Diff+cyt+slow 2", "Sim: Diff+cyt+slow 3", "Synth: ER", "Synth: SF 1", "Synth: SF 2", "Synth: SF 3", "Synth: WS 1", "Synth: WS 2", "Synth: WS 3", "Synth: Geometric", "Synth: Star")
for(cs in (3+seq(from=0,to=7)*5)) {
  for(ct in 1:2) {
      all.titles = c(all.titles, paste("Clique: ", cs, "-", ct, sep=""))
    }
}
cbPalette = c(rep("#FFAAAA", 16), rep("#AAAAFF", 12), rep("#AAFFAA", 9), "#000000")

my.labels = sort(all.titles)
my.labels[length(my.labels)] = "Bio"

subset$L = factor(subset$L)
subset$expt = all.titles[subset$expt+1]
png(output.fname, width=1200,height=1200)
ggplot(data = subset, aes(x=prop.edges,y=prop.bingos,group=expt)) +
  stat_smooth(aes(group=expt, color= expt, fill=expt), method="loess") +
  facet_wrap(~L, scales="free") +
  scale_fill_manual(labels = my.labels, values=cbPalette) +
  scale_colour_manual(labels = my.labels, values=cbPalette) +
#  scale_fill_manual(values=cbPalette) +
#  scale_colour_manual(values=cbPalette) +
  xlab("Proportion of edges used") + ylab("Bingo score") + labs(colour = "Network", fill="Network") +
  theme(axis.title.x = element_text(size=24), axis.title.y = element_text(size=24),
        axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
	legend.title = element_text(size=24), legend.text = element_text(size=14),
	strip.text.x = element_text(size = 18))
  
dev.off()

