source("/Users/ogunadebali/SancarLab/scripts/tcsv.R")

mappedReads <- function(file){
  mappedReadsDf <- read.tcsv(file, sep="\t")
  return(mappedReadsDf)
}

