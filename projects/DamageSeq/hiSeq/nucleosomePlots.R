library("ggplot2")
source("/Users/ogunadebali/SancarLab/scripts/tcsv.R")
source("/Users/ogunadebali/SancarLab/scripts/multiplot.R")
source("/Users/ogunadebali/SancarLab/scripts/mappedReads.R")

  keyWord = "Data"
  keyWord = "DataCopy"

  keyWord = "DataTop"
  keyWord = "DataLow"
  keyWord = "DataRep1T1M"
  keyWord = "DataRep1copy"
  keyWord = "Data_rep2"
  input = paste("/Users/ogunadebali/SancarLab/DamageSeq/hiSeq/mergedNucleosome", keyWord, ".csv", sep="")
  input = paste("/Users/ogunadebali/SancarLab/DamageSeq/hiSeq/mergedNucleosomeData_rep.intersect.csv", sep="")
  input = paste("/Users/ogunadebali/SancarLab/DamageSeq/hiSeq/mergedNucleosomeData_rep.intersect_rep3.csv", sep="")
  input = paste("/Users/ogunadebali/SancarLab/DamageSeq/hiSeq/mergedNucleosomeData_rep.intersect_adar.csv", sep="")
  input = paste("/Users/ogunadebali/SancarLab/DamageSeq/hiSeq/mergedNuPeaks_NHF1.csv", sep="")
  

  df <- read.tcsv(input)
  
  df$positions <- c(-73:73)
  #df$positions <- c(-70:70)
  
  
  normalize <- function(data) {
    result <- (data / sum(data)) * 1000000
    return(result)
  }
  df$GM12878_6.4_20J_cell_A_Pls <- normalize(df$GM12878_6.4_20J_cell_A_Pls)
  df$GM12878_6.4_20J_cell_A_Min <- normalize(df$GM12878_6.4_20J_cell_A_Min)
  df$GM12878_6.4_20J_cell_B_Pls <- normalize(df$GM12878_6.4_20J_cell_B_Pls)
  df$GM12878_6.4_20J_cell_B_Min <- normalize(df$GM12878_6.4_20J_cell_B_Min)

  df$GM12878_6.4_20J_nakedDNA_A_Pls <- normalize(df$GM12878_6.4_20J_nakedDNA_A_Pls)
  df$GM12878_6.4_20J_nakedDNA_A_Min <- normalize(df$GM12878_6.4_20J_nakedDNA_A_Min)
  df$GM12878_6.4_20J_nakedDNA_B_Pls <- normalize(df$GM12878_6.4_20J_nakedDNA_B_Pls)
  df$GM12878_6.4_20J_nakedDNA_B_Min <- normalize(df$GM12878_6.4_20J_nakedDNA_B_Min)
  
  df$GM12878_CPD_20J_cell_A_Pls <- normalize(df$GM12878_CPD_20J_cell_A_Pls)
  df$GM12878_CPD_20J_cell_A_Min <- normalize(df$GM12878_CPD_20J_cell_A_Min)
  df$GM12878_CPD_20J_cell_B_Pls <- normalize(df$GM12878_CPD_20J_cell_B_Pls)
  df$GM12878_CPD_20J_cell_B_Min <- normalize(df$GM12878_CPD_20J_cell_B_Min)
  
  df$GM12878_CPD_20J_nakedDNA_A_Pls <- normalize(df$GM12878_CPD_20J_nakedDNA_A_Pls)
  df$GM12878_CPD_20J_nakedDNA_A_Min <- normalize(df$GM12878_CPD_20J_nakedDNA_A_Min)
  df$GM12878_CPD_20J_nakedDNA_B_Pls <- normalize(df$GM12878_CPD_20J_nakedDNA_B_Pls)
  df$GM12878_CPD_20J_nakedDNA_B_Min <- normalize(df$GM12878_CPD_20J_nakedDNA_B_Min)
  
  
  df$GM12878_Cisplatin_Cell_A_Pls <- normalize(df$GM12878_Cisplatin_Cell_A_Pls)
  df$GM12878_Cisplatin_Cell_A_Min <- normalize(df$GM12878_Cisplatin_Cell_A_Min)
  df$GM12878_Cisplatin_Cell_B_Pls <- normalize(df$GM12878_Cisplatin_Cell_B_Pls)
  df$GM12878_Cisplatin_Cell_B_Min <- normalize(df$GM12878_Cisplatin_Cell_B_Min)
  
  df$GM12878_Cisplatin_nakedDNA_A_Pls <- normalize(df$GM12878_Cisplatin_nakedDNA_A_Pls)
  df$GM12878_Cisplatin_nakedDNA_A_Min <- normalize(df$GM12878_Cisplatin_nakedDNA_A_Min)
  df$GM12878_Cisplatin_nakedDNA_B_Pls <- normalize(df$GM12878_Cisplatin_nakedDNA_B_Pls)
  df$GM12878_Cisplatin_nakedDNA_B_Min <- normalize(df$GM12878_Cisplatin_nakedDNA_B_Min)
  
  
  
  
  #df$GM12878_6.4_20J_cell_A_Pls
  #df$GM12878_6.4_20J_cell_B_Pls
  df$GM12878_6.4_20J_cell_A_Min_reversed <- rev(df$GM12878_6.4_20J_cell_A_Min);
  df$GM12878_6.4_20J_cell_B_Min_reversed <- rev(df$GM12878_6.4_20J_cell_B_Min);
  
  #df$GM12878_6.4_20J_nakedDNA_A_Pls
  df$GM12878_6.4_20J_nakedDNA_A_Min_reversed <- rev(df$GM12878_6.4_20J_nakedDNA_A_Min);
  #df$GM12878_6.4_20J_nakedDNA_B_Pls
  df$GM12878_6.4_20J_nakedDNA_B_Min_reversed <- rev(df$GM12878_6.4_20J_nakedDNA_B_Min);
  
  #df$GM12878_CPD_20J_cell_A_Pls
  df$GM12878_CPD_20J_cell_A_Min_reversed <- rev(df$GM12878_CPD_20J_cell_A_Min);
  #df$GM12878_CPD_20J_cell_B_Pls
  df$GM12878_CPD_20J_cell_B_Min_reversed <- rev(df$GM12878_CPD_20J_cell_B_Min);
  
  #df$GM12878_CPD_20J_nakedDNA_A_Pls
  df$GM12878_CPD_20J_nakedDNA_A_Min_reversed <- rev(df$GM12878_CPD_20J_nakedDNA_A_Min);
  #df$GM12878_CPD_20J_nakedDNA_B_Pls
  df$GM12878_CPD_20J_nakedDNA_B_Min_reversed <- rev(df$GM12878_CPD_20J_nakedDNA_B_Min);
  
  df$GM12878_Cisplatin_cell_A_Min_reversed <- rev(df$GM12878_Cisplatin_Cell_A_Min)
  df$GM12878_Cisplatin_nakedDNA_A_Min_reversed <- rev(df$GM12878_Cisplatin_nakedDNA_A_Min)
  df$GM12878_Cisplatin_cell_B_Min_reversed <- rev(df$GM12878_Cisplatin_Cell_B_Min)
  df$GM12878_Cisplatin_nakedDNA_B_Min_reversed <- rev(df$GM12878_Cisplatin_nakedDNA_B_Min)
  
  sets = c(
    c("GM12878_CPD_20J_cell_A_Min_reversed", "GM12878_CPD_20J_nakedDNA_A_Min_reversed", "GM12878_CPD_20J_cell_A_Pls", "GM12878_CPD_20J_nakedDNA_A_Pls"),
    c("GM12878_CPD_20J_cell_B_Min_reversed", "GM12878_CPD_20J_nakedDNA_B_Min_reversed", "GM12878_CPD_20J_cell_B_Pls", "GM12878_CPD_20J_nakedDNA_B_Pls"),
    c("GM12878_6.4_20J_cell_A_Min_reversed", "GM12878_6.4_20J_nakedDNA_A_Min_reversed", "GM12878_6.4_20J_cell_A_Pls", "GM12878_6.4_20J_nakedDNA_A_Pls"),
    c("GM12878_6.4_20J_cell_B_Min_reversed", "GM12878_6.4_20J_nakedDNA_B_Min_reversed", "GM12878_6.4_20J_cell_B_Pls", "GM12878_6.4_20J_nakedDNA_B_Pls"),
    c("GM12878_Cisplatin_cell_A_Min_reversed", "GM12878_Cisplatin_nakedDNA_A_Min_reversed", "GM12878_Cisplatin_Cell_A_Pls", "GM12878_Cisplatin_nakedDNA_A_Pls"),
    c("GM12878_Cisplatin_cell_B_Min_reversed", "GM12878_Cisplatin_nakedDNA_B_Min_reversed", "GM12878_Cisplatin_Cell_B_Pls", "GM12878_Cisplatin_nakedDNA_B_Pls")
    )
  
  headers = c("CPD_A", "CPD_B", "6-4_A", "6-4_B", "Cisplatin_A", "Cisplatin_B" )
  
  #linePlot <- function(i, ymin=0.9, ymax=1.3) {
linePlot <- function(i, ymin=0.9, ymax=1.3) {
      
    p <- ggplot(df) + 
    #geom_line(aes(x=positions, y= df[sets[i*4 + 1]]), colour= "brown", size=1.1) +
      geom_line(aes(x=positions, y= df[sets[i*4 + 1]]/df[sets[i*4 + 2]]), colour= "brown", size=1.1) +
      geom_line(aes(x=positions, y=df[sets[i*4 + 3]]/df[sets[i*4 + 4]]), colour= "green", size=1.1) +
    #geom_line(aes(x=positions, y=df[sets[i*4 + 2]]), colour= "green", size=1.1) +
      
    #ylim(ymin, ymax) +
    #ggtitle(headers[i+1]) +
    ylab("Cell / Naked Damage") +
    xlab("Position Relative to Nucleosome Dyad") + 
    labs(title = headers[i+1]) +
    theme(plot.title = element_text(hjust = 0.5))
    return(p)
  }
  
  
  #pdf(paste("~/SancarLab/DamageSeq/hiSeq/nucleosome_CPD_", keyWord, ".pdf", sep=""))
  multiplot(linePlot(0, 0.90, 1.06), linePlot(1, 0.92, 1.06), cols = 1)
  #dev.off()
  #pdf(paste("~/SancarLab/DamageSeq/hiSeq/nucleosome_64_", keyWord, ".pdf", sep=""))
  multiplot(linePlot(2, 0.97, 1.05), linePlot(3,  0.89, 1.07), cols = 1)
  #  dev.off()
  #pdf(paste("~/SancarLab/DamageSeq/hiSeq/nucleosome_Cisplatin_", keyWord, ".pdf", sep=""))
  multiplot(linePlot(4, 0.97, 1.1), linePlot(5, 0.97, 1.25), cols = 1)
  #dev.off()


multiplot(linePlot(0), linePlot(1), linePlot(2), linePlot(3), linePlot(4), linePlot(5), cols=2)

linePlot <- function(i) {
  p <- ggplot(df) + 
    geom_line(aes(x=positions, y= df[sets[i + 1]]), colour= "brown", size=1.1) +
    geom_line(aes(x=positions, y= df[sets[i + 3]]), colour= "green", size=1.1) +
    ylab("Damage Count") +
    xlab("Position Relative to Nucleosome Dyad") + 
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

multiplot(linePlot(0), linePlot(4), linePlot(1), linePlot(5), cols=2)
multiplot(linePlot(8), linePlot(12), linePlot(9), linePlot(13), cols=2)
multiplot(linePlot(16), linePlot(20), linePlot(17), linePlot(21), cols=2)

