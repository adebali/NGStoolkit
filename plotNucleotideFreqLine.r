#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
input = args[1]
treatment = args[2]
output = paste(input,'.',treatment,'.pdf',sep='')


pdf(output,width=10,height=10)
dat <- as.matrix(read.table(input, header = TRUE))
title = gsub("_", " ", paste(treatment, "-", "Dinucleotide Frequencies"))
positionLength = length(dat[1,])
colors = rev(c("#a1dab4", 
               "#41b6c4", 
               "#2c7fb8", 
               "#253494"))
colors16 = c("#a6cee3",
              "#1f78b4",
              "#b2df8a",
              "#33a02c",
              "#fb9a99",
              "#e31a1c",
              "#fdbf6f",
              "#ff7f00",
              "#cab2d6",
              "#6a3d9a",
              "#ffff99",
              "#b15928",
              "#555577",
              "#036d24",
              "#cb2f2a",
              "#972acb")
#library(RColorBrewer)
#colors = rev(brewer.pal(6,"Blues"))[1:4]


# Get frequencies
for (i in 1:positionLength){
  dat[,i] = dat[,i]/sum(dat[,i])
}


par(fig=c(0.1,0.95,0.1,0.95))
plot(NA, xlim=c(1,positionLength), NA,  ylim=c(0,1), main=title, bty='n', xlab='Position', ylab='Frequency')

for (i in 1:length(dat[,1])){
  points(1:positionLength, as.list(dat[i,]), 
         type="b", 
         lty=i, 
         pch=i,
         col=colors16[i],
         bg = colors16[i])
}

par(fig=c(0.7,0.9,0.1,0.95), new=TRUE)

legend( x="topright", 
        legend=row.names(dat),
        col=colors16, lwd=1, lty=c(1,16), 
        pch=c(1:16), bty = 'n', cex= 1)

dev.off()