#!/usr/bin/env Rscript
library(dplyr)
library(reshape2)

setwd("/proj/sancarlb/users/ogun/sancarlabutils/projects/replication")
input = 'dataDir/merged_1kCounts.txt'
dat <- read.csv(input, header = T, sep="\t")
filtered <- filter(dat, chromosome == "Chr2", start > 113000000, start < 130000000)
filtered_dcasted <- dcast(filtered, start + product + replicate + phase + method ~ strandName, value.var="count")
filtered_dcasted$direction <- (filtered_dcasted$plus - filtered_dcasted$minus) / (filtered_dcasted$plus + filtered_dcasted$minus)
write.csv(filtered_dcasted, "filtered.csv", row.names = F)
