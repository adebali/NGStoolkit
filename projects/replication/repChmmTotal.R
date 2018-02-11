args = commandArgs(trailingOnly=TRUE)
d <- read.csv(args[1], header = F, sep="\t")
names(d) <- c("chromosome", "start", "end", "name", "score", "strand", "count")
d$length <- d$end - d$start
total <- sum(d$count)

library(plyr)
dagg <- ddply(d, .(name), 
      summarize, 
      count = sum(count),
      length = sum(length))
dagg$RPM <- ((dagg$count/total)*1e6)
dagg$RPKM <- (dagg$RPM/dagg$length)*1e3

getField <- function(t, field) {
  return (strsplit(t, '_')[[1]][field])
}

dagg$RD <- lapply(as.character(dagg$name), getField, 1)
dagg$CS <- lapply(as.character(dagg$name), getField, 2)
library(dplyr)
dat <- select(dagg, RD, CS, length, RPM, RPKM)
df <- apply(dat,2,as.character)
write.table(df, file = args[2], sep = '\t', quote = F, row.names = F, col.names = F)
