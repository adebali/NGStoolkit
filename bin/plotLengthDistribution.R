#!/usr/bin/env Rscript

library("ggplot2")
library("reshape2")
source("~/ogun/scripts/tcsv.R")

args<-commandArgs(TRUE)
input = args[1]
treatment = args[2]
output = paste(input,'_lenDist','.pdf',sep='')
dat <- read.tcsv(input, header = T,check.names=F, sep="\t")

pdf(output)
ggplot(melt(dat)) + geom_bar(stat="identity", aes(x=variable, y=value)) + ggtitle(treatment) +
ylab("Read count") + xlab("Position")
dev.off()