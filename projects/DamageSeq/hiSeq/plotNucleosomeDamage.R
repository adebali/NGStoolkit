#!/usr/bin/env Rscript

library("ggplot2")
data = read.table(file("stdin"), sep="\t", header = F)
args<-commandArgs(TRUE)
output = args[1]
multiplicationFactor = as.numeric(args[2])
treatment = args[3]
oneSideFlanking = as.numeric(args[4])
counts <- table(data) * multiplicationFactor
negativeSide = -1 * oneSideFlanking
positiveSide = oneSideFlanking
positions <- c(negativeSide:positiveSide)
d = data.frame(positions, counts)
print(d)
lapply(counts, write, output, append=T, ncolumns=1)

# ## Draw an individual plot
# pdf(paste(output, ".pdf", sep=""),width=5,height=5)
# p <- ggplot(d, aes(x=positions, y=Freq))
# p + geom_line() + ggtitle(treatment)
# dev.off()


