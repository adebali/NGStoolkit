#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
input = args[1]
treatment = args[2]
output = paste(input,'.',treatment,'.chromHMM.pdf',sep='')
pdf(output,width=5,height=5)
title = gsub("_", " ", paste(treatment, "-", "ChromHMM"))

data = bedCountTable = read.table(input)
colors =  c("red","pink","purple",
            "orange","orange","yellow","yellow",
            "royalblue4",
            "chartreuse4","chartreuse4", "olivedrab1", 
            "gray50", "gray50",
            "gray75", "gray75")

labels = c("Active promoter", "Waek promoter", "Poised promoter",
           "Strong enhancer", "Strong enhancer", "Weak enhancer", "Weak enhancer",
           "Insulator",
           "Txn transition", "Txn elongation", "Weak txn",
           "Repressed", "Heterochrom",
           "Repetetive", "Repetetive"
           )

par(fig=c(0.05,0.95,0.2,0.95))
boxplot(data, outline=FALSE, 
      col=(colors), cex=0.1, 
      names=labels, las=2, 
      ylim=c(0,30),
      main=title, border= rgb(0.2,0.2,0.2,1), frame=F)

dev.off()
