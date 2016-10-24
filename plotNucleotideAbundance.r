#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
input = args[1]
treatment = args[2]
output = paste(input,'.',treatment,'.pdf',sep='')

pdf(output,width=5,height=5)
dat <- as.matrix(read.table(input, header = TRUE))
title = gsub("_", " ", paste(treatment, "-", "Nucleotide Frequencies"))
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
for (i in 1:length(colnames(dat))){
    dat[,i] = dat[,i]/sum(dat[,i])
  }

par(xpd=T, mar=par()$mar+c(0,0,0,4));

barplot(dat, 
        border= FALSE,
        space = c(0.4,0.4),
        main = title, 
        xlab="Position", 
        ylab="Frequency",
        names.arg=c(1:positionLength),
        col=colors,
        cex.names=1);

legend(positionLength + 5, 0.75,
       rev(row.names(dat)), 
       cex=1, 
       fill=rev(colors), 
       border = FALSE,
       y.intersp = 2,
       bty='n');
  
# Restore default clipping rect
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()