#!/usr/bin/env Rscript
lengths <- as.integer(readLines("stdin"))
args<-commandArgs(TRUE)
output = args[1]
output = paste(output,'.pdf',sep='')
pdf(output,width=5,height=5)
hist(lengths, seq(0,10000,by=100))
dev.off()