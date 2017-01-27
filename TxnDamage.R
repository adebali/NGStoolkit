library("ggplot2")
source("/Users/ogunadebali/SancarLab/scripts/tcsv.R")
source("/Users/ogunadebali/SancarLab/scripts/multiplot.R")


#input = paste("/Users/ogunadebali/SancarLab/DamageSeq/hiSeq/mergedTxn.csv", sep="")
input = paste("/Users/ogunadebali/SancarLab/DamageSeq/hiSeq/mergedTFBS.csv", sep="")

df <- read.tcsv(input)

df$positions <- c(-1000:1000)

## 6-4
# Cell
df$GM12878_6.4_20J_cell_A_Pls <- (df$GM12878_6.4_20J_cell_A_Pls / sum(df$GM12878_6.4_20J_cell_A_Pls)) * 1000000
df$GM12878_6.4_20J_cell_B_Pls <- (df$GM12878_6.4_20J_cell_B_Pls / sum(df$GM12878_6.4_20J_cell_B_Pls)) * 1000000
df$GM12878_6.4_20J_cell_A_Min <- (df$GM12878_6.4_20J_cell_A_Min / sum(df$GM12878_6.4_20J_cell_A_Min)) * 1000000
df$GM12878_6.4_20J_cell_B_Min <- (df$GM12878_6.4_20J_cell_B_Min / sum(df$GM12878_6.4_20J_cell_B_Min)) * 1000000
df$GM12878_6.4_20J_cell_A_Min_reversed <- rev(df$GM12878_6.4_20J_cell_A_Min);
df$GM12878_6.4_20J_cell_B_Min_reversed <- rev(df$GM12878_6.4_20J_cell_B_Min);

# Naked
df$GM12878_6.4_20J_nakedDNA_A_Pls <- (df$GM12878_6.4_20J_nakedDNA_A_Pls / sum(df$GM12878_6.4_20J_nakedDNA_A_Pls)) * 1000000
df$GM12878_6.4_20J_nakedDNA_B_Pls <- (df$GM12878_6.4_20J_nakedDNA_B_Pls / sum(df$GM12878_6.4_20J_nakedDNA_B_Pls)) * 1000000
df$GM12878_6.4_20J_nakedDNA_A_Min <- (df$GM12878_6.4_20J_nakedDNA_A_Min / sum(df$GM12878_6.4_20J_nakedDNA_A_Min)) * 1000000
df$GM12878_6.4_20J_nakedDNA_B_Min <- (df$GM12878_6.4_20J_nakedDNA_B_Min / sum(df$GM12878_6.4_20J_nakedDNA_B_Min)) * 1000000
df$GM12878_6.4_20J_nakedDNA_A_Min_reversed <- rev(df$GM12878_6.4_20J_nakedDNA_A_Min);
df$GM12878_6.4_20J_nakedDNA_B_Min_reversed <- rev(df$GM12878_6.4_20J_nakedDNA_B_Min);
##

## CPD
# Cell
df$GM12878_CPD_20J_cell_A_Pls <- (df$GM12878_CPD_20J_cell_A_Pls / sum(df$GM12878_CPD_20J_cell_A_Pls)) * 1000000
df$GM12878_CPD_20J_cell_B_Pls <- (df$GM12878_CPD_20J_cell_B_Pls / sum(df$GM12878_CPD_20J_cell_B_Pls)) * 1000000
df$GM12878_CPD_20J_cell_A_Min <- (df$GM12878_CPD_20J_cell_A_Min / sum(df$GM12878_CPD_20J_cell_A_Min)) * 1000000
df$GM12878_CPD_20J_cell_B_Min <- (df$GM12878_CPD_20J_cell_B_Min / sum(df$GM12878_CPD_20J_cell_B_Min)) * 1000000
df$GM12878_CPD_20J_cell_A_Min_reversed <- rev(df$GM12878_CPD_20J_cell_A_Min);
df$GM12878_CPD_20J_cell_B_Min_reversed <- rev(df$GM12878_CPD_20J_cell_B_Min);
# Naked
df$GM12878_CPD_20J_nakedDNA_A_Pls <- (df$GM12878_CPD_20J_nakedDNA_A_Pls / sum(df$GM12878_CPD_20J_nakedDNA_A_Pls)) * 1000000
df$GM12878_CPD_20J_nakedDNA_B_Pls <- (df$GM12878_CPD_20J_nakedDNA_B_Pls / sum(df$GM12878_CPD_20J_nakedDNA_B_Pls)) * 1000000
df$GM12878_CPD_20J_nakedDNA_A_Min <- (df$GM12878_CPD_20J_nakedDNA_A_Min / sum(df$GM12878_CPD_20J_nakedDNA_A_Min)) * 1000000
df$GM12878_CPD_20J_nakedDNA_B_Min <- (df$GM12878_CPD_20J_nakedDNA_B_Min / sum(df$GM12878_CPD_20J_nakedDNA_B_Min)) * 1000000
df$GM12878_CPD_20J_nakedDNA_A_Min_reversed <- rev(df$GM12878_CPD_20J_nakedDNA_A_Min);
df$GM12878_CPD_20J_nakedDNA_B_Min_reversed <- rev(df$GM12878_CPD_20J_nakedDNA_B_Min);

## Cisplatin
# Cell
df$GM12878_Cisplatin_Cell_A_Pls <- (df$GM12878_Cisplatin_Cell_A_Pls / sum(df$GM12878_Cisplatin_Cell_A_Pls)) * 1000000
df$GM12878_Cisplatin_Cell_B_Pls <- (df$GM12878_Cisplatin_Cell_B_Pls / sum(df$GM12878_Cisplatin_Cell_B_Pls)) * 1000000
df$GM12878_Cisplatin_Cell_A_Min <- (df$GM12878_Cisplatin_Cell_A_Min / sum(df$GM12878_Cisplatin_Cell_A_Min)) * 1000000
df$GM12878_Cisplatin_Cell_B_Min <- (df$GM12878_Cisplatin_Cell_B_Min / sum(df$GM12878_Cisplatin_Cell_B_Min)) * 1000000
df$GM12878_Cisplatin_Cell_A_Min_reversed <- rev(df$GM12878_Cisplatin_Cell_A_Min);
df$GM12878_Cisplatin_Cell_B_Min_reversed <- rev(df$GM12878_Cisplatin_Cell_B_Min);
# Naked
df$GM12878_Cisplatin_nakedDNA_A_Pls <- (df$GM12878_Cisplatin_nakedDNA_A_Pls / sum(df$GM12878_Cisplatin_nakedDNA_A_Pls)) * 1000000
df$GM12878_Cisplatin_nakedDNA_B_Pls <- (df$GM12878_Cisplatin_nakedDNA_B_Pls / sum(df$GM12878_Cisplatin_nakedDNA_B_Pls)) * 1000000
df$GM12878_Cisplatin_nakedDNA_A_Min <- (df$GM12878_Cisplatin_nakedDNA_A_Min / sum(df$GM12878_Cisplatin_nakedDNA_A_Min)) * 1000000
df$GM12878_Cisplatin_nakedDNA_B_Min <- (df$GM12878_Cisplatin_nakedDNA_B_Min / sum(df$GM12878_Cisplatin_nakedDNA_B_Min)) * 1000000
df$GM12878_Cisplatin_nakedDNA_A_Min_reversed <- rev(df$GM12878_Cisplatin_nakedDNA_A_Min);
df$GM12878_Cisplatin_nakedDNA_B_Min_reversed <- rev(df$GM12878_Cisplatin_nakedDNA_B_Min);


sets = c(
  c("GM12878_CPD_20J_cell_A_Min_reversed", "GM12878_CPD_20J_nakedDNA_A_Min_reversed", "GM12878_CPD_20J_cell_A_Pls", "GM12878_CPD_20J_nakedDNA_A_Pls"),
  c("GM12878_CPD_20J_cell_B_Min_reversed", "GM12878_CPD_20J_nakedDNA_B_Min_reversed", "GM12878_CPD_20J_cell_B_Pls", "GM12878_CPD_20J_nakedDNA_B_Pls"),
  c("GM12878_6.4_20J_cell_A_Min_reversed", "GM12878_6.4_20J_nakedDNA_A_Min_reversed", "GM12878_6.4_20J_cell_A_Pls", "GM12878_6.4_20J_nakedDNA_A_Pls"),
  c("GM12878_6.4_20J_cell_B_Min_reversed", "GM12878_6.4_20J_nakedDNA_B_Min_reversed", "GM12878_6.4_20J_cell_B_Pls", "GM12878_6.4_20J_nakedDNA_B_Pls"),
  c("GM12878_Cisplatin_Cell_A_Min_reversed", "GM12878_Cisplatin_nakedDNA_A_Min_reversed", "GM12878_Cisplatin_Cell_A_Pls", "GM12878_Cisplatin_nakedDNA_A_Pls"),
  c("GM12878_Cisplatin_Cell_B_Min_reversed", "GM12878_Cisplatin_nakedDNA_B_Min_reversed", "GM12878_Cisplatin_Cell_B_Pls", "GM12878_Cisplatin_nakedDNA_B_Pls")
)

headers = c("CPD_A", "CPD_B", "6-4_A", "6-4_B", "Cisplatin_A", "Cisplatin_B" )

linePlot <- function(i, ymin=0.9, ymax=1.3) {
  p <- ggplot(df) + 
    geom_line(aes(x=positions, y= df[sets[i*4 + 1]]/df[sets[i*4 + 2]]), colour= "brown", size=0.2) +
    geom_line(aes(x=positions, y=df[sets[i*4 + 3]]/df[sets[i*4 + 4]]), colour= "green", size=0.2) +
    #ylim(ymin, ymax) +
    #ggtitle(headers[i+1]) +
    ylab("Cell / Naked Damage") +
    xlab("Position Relative to TBS Center") + 
    labs(title = headers[i+1]) +
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

#pdf(paste("~/SancarLab/DamageSeq/hiSeq/tTxn_CPD_", keyWord, ".pdf", sep=""))
multiplot(linePlot(0, 0.90, 1.06), linePlot(1, 0.92, 1.06), cols = 1)
#dev.off()
#pdf(paste("~/SancarLab/DamageSeq/hiSeq/Txn_64_", keyWord, ".pdf", sep=""))
multiplot(linePlot(2, 0.97, 1.05), linePlot(3,  0.89, 1.07), cols = 1)
#dev.off()
#pdf(paste("~/SancarLab/DamageSeq/hiSeq/Txn_Cisplatin_", keyWord, ".pdf", sep=""))
multiplot(linePlot(4, 0.97, 1.1), linePlot(5, 0.97, 1.25), cols = 1)
#dev.off()

