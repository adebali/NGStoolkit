#!/usr/bin/env Rscript
library("ggplot2")
library("reshape2")


input = "dataDir/merged_geneCounts.txt"

d <- read.csv(input, header = T, sep="\t")

d$treatment_title <- factor(d$treatment_title, levels=unique(d[order(d$no),]$treatment_title))
d$ZT <- factor(d$ZT, levels=unique(d[order(d$no),]$ZT))
d$organ <- factor(d$organ, levels=unique(d[order(d$no),]$organ))
d$TSNTS <- factor(d$TSNTS, levels=c("TS", "NTS"))

# pdf("geneCount_points.pdf")
# ggplot(d, aes(x=ZT, y=count), log10="y") + geom_jitter(width = 0.25, aes(color=TSNTS), size=0.1, alpha= 0.01) + facet_grid(~organ~TSNTS) + scale_y_log10()
# ggplot(d, aes(x=ZT, y=count), log10="y") + geom_jitter(width = 0.25, aes(color=TSNTS), size=0.1, alpha= 0.01) +  geom_boxplot(outlier.colour=NA, fill=NA, colour="grey20") + facet_grid(~organ~TSNTS) + scale_y_log10()
# dev.off()

h <- head(d, n =605713)
mean <- mean(h$count)
variance <- var(h$count)
h$zscore <- (h$count - mean)/variance

# d$zscore <- (d$count - mean)/variance
# d$no <- rep(list(c(1:length(d[,1])/24)), 24)

head(h)
# pdf("heatmap.pdf")
# ggplot(h, aes(ZT, no)) + geom_tile(aes(fill=zscore)) + scale_fill_gradient(low = "blue", high = "yellow") + facet_grid(~organ~TSNTS) 
# dev.off()

# d$ZT_ = paste("ZT", d$ZT, sep = "")

# t <- dcast(d[which(d$TSNTS=="TS" & d$organ=="liver"),], no + chr + start + end + name ~ ZT_, value.var="count")
# tnz <- t[which(t$ZT0 != 0 & t$ZT0 != 0 & t$ZT4 != 0 & t$ZT8 != 0 & t$ZT12 != 0 & t$ZT16 != 0 & t$ZT20 != 0),]
# tnzm <- melt(tnz, id=c("no","chr","start","end"))


