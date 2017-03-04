library("ggplot2")
library("plotly")
source("/Users/ogunadebali/SancarLab/scripts/tcsv.R")
source("/Users/ogunadebali/SancarLab/scripts/multiplot.R")


#input = paste("/Users/ogunadebali/SancarLab/DamageSeq/hiSeq/mergedTxn.csv", sep="")
input = paste("/Users/ogunadebali/SancarLab/DamageSeq/hiSeq/mergedTFBS.csv", sep="")
input = paste("/Users/ogunadebali/SancarLab/DamageSeq/hiSeq/mergedTFBS_GM12878.csv", sep="")

nucInput <- "/Users/ogunadebali/SancarLab/DamageSeq/hiSeq/Dinucleotides/wgEncodeRegTfbsClusteredWithCellsV3.GM12878.sorted.1K.sorted.fa.freq.csv"
freq <- read.tcsv(nucInput, header = T, sep="\t")
freq$positions <- c(-999:999)
df <- read.tcsv(input)

df$AA<-c(NaN, freq$AA, NaN)
df$TT<-c(NaN, freq$TT, NaN)
df$GG<-c(NaN, freq$GG, NaN)
df$CC<-c(NaN, freq$CC, NaN)
df$AT<-c(NaN, freq$AT, NaN)
df$TA<-c(NaN, freq$TA, NaN)
df$AG<-c(NaN, freq$AG, NaN)
df$GA<-c(NaN, freq$GA, NaN)
df$AC<-c(NaN, freq$AC, NaN)
df$CA<-c(NaN, freq$CA, NaN)
df$GC<-c(NaN, freq$GC, NaN)
df$CG<-c(NaN, freq$CG, NaN)
df$TC<-c(NaN, freq$TC, NaN)
df$CT<-c(NaN, freq$CT, NaN)
df$GT<-c(NaN, freq$GT, NaN)
df$TG<-c(NaN, freq$TG, NaN)

df$positions <- c(-1000:1000)


#p1 <- ggplot() + 
#  geom_line(data=df, aes(x=positions, y= GM12878_6.4_20J_cell_A_Pls), colour= "brown", size=0.2)

#p2 <- ggplot() + 
#  geom_line(data=df, aes(x=positions, y= AA), colour= "purple", size=0.2)


#library(grid)
#grid.newpage()

#grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))

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
df$GM12878_6.4_20J_nakedDNA_B_Pls <- (df$GM12878_6.4_20J_nakedDNA_C_Pls / sum(df$GM12878_6.4_20J_nakedDNA_C_Pls)) * 1000000
df$GM12878_6.4_20J_nakedDNA_A_Min <- (df$GM12878_6.4_20J_nakedDNA_A_Min / sum(df$GM12878_6.4_20J_nakedDNA_A_Min)) * 1000000
df$GM12878_6.4_20J_nakedDNA_B_Min <- (df$GM12878_6.4_20J_nakedDNA_C_Min / sum(df$GM12878_6.4_20J_nakedDNA_C_Min)) * 1000000
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
df$GM12878_CPD_20J_nakedDNA_B_Pls <- (df$GM12878_CPD_20J_nakedDNA_C_Pls / sum(df$GM12878_CPD_20J_nakedDNA_C_Pls)) * 1000000
df$GM12878_CPD_20J_nakedDNA_A_Min <- (df$GM12878_CPD_20J_nakedDNA_A_Min / sum(df$GM12878_CPD_20J_nakedDNA_A_Min)) * 1000000
df$GM12878_CPD_20J_nakedDNA_B_Min <- (df$GM12878_CPD_20J_nakedDNA_C_Min / sum(df$GM12878_CPD_20J_nakedDNA_C_Min)) * 1000000
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


setsHeaders = c(
  c("CPD_cell_A_Min", "CPD_nakedDNA_A_Min", "CPD_cell_A_Pls", "CPD_nakedDNA_A_Pls"),
  c("CPD_cell_B_Min", "CPD_nakedDNA_B_Min", "CPD_cell_B_Pls", "CPD_nakedDNA_B_Pls"),
  c("6.4_cell_A_Min", "6.4_nakedDNA_A_Min", "6.4_cell_A_Pls", "6.4_nakedDNA_A_Pls"),
  c("6.4_cell_B_Min", "6.4_nakedDNA_B_Min", "6.4_cell_B_Pls", "6.4_nakedDNA_B_Pls"),
  c("Cisplatin_Cell_A_Min", "Cisplatin_nakedDNA_A_Min", "Cisplatin_Cell_A_Pls", "Cisplatin_nakedDNA_A_Pls"),
  c("Cisplatin_Cell_B_Min", "Cisplatin_nakedDNA_B_Min", "Cisplatin_Cell_B_Pls", "Cisplatin_nakedDNA_B_Pls")
)


headers = c("CPD_A", "CPD_B", "6-4_A", "6-4_B", "Cisplatin_A", "Cisplatin_B" )

linePlot <- function(i, ymin=0.9, ymax=1.3) {
  p <- ggplot(df) + 
    geom_smooth(aes(x=positions, y= df[sets[i*4 + 1]]/df[sets[i*4 + 2]]), colour= "brown", size=0.2) +
    geom_smooth(aes(x=positions, y=df[sets[i*4 + 3]]/df[sets[i*4 + 4]]), colour= "green", size=0.2) +
    ylim(ymin, ymax) +
    ylab("Cell / Naked Damage") +
    xlab("Position Relative to TBS Center") + 
    labs(title = headers[i+1]) +
    theme(plot.title = element_text(hjust = 0.5))
  # return(p)
  ggplotly(p)
}

linePlot(5)

#pdf(paste("~/SancarLab/DamageSeq/hiSeq/TFBS_CPD.pdf", sep=""))
multiplot(linePlot(0,  0.85, 1.1), linePlot(1, 0.85, 1.1), cols = 1)
#dev.off()
#pdf(paste("~/SancarLab/DamageSeq/hiSeq/TFBS_64.pdf", sep=""))
multiplot(linePlot(2, 0.85, 1.1), linePlot(3, 0.85, 1.1), cols = 1)
#dev.off()
#pdf(paste("~/SancarLab/DamageSeq/hiSeq/TFBS_Cisplatin.pdf", sep=""))
multiplot(linePlot(4, 0.85, 1.1), linePlot(5, 0.85, 1.1), cols = 1)
#dev.off()

linePlotUV <- function(i) {
  p <- ggplot(df) + 
    geom_line(aes(x=positions, y= df[sets[i]]/ (df$TT + df$AA + df$TC + df$GA + )), colour= "brown", size=0.1) +
    #geom_line(aes(x=positions, y=df[sets[i*4 + 3]]/df[sets[i*4 + 4]]), colour= "green", size=0.2) +
    #ylim(350, 600) +
    ylab("Damage Per 1M Mapped Reads") +
    xlab("Position Relative to TBS Center") +
    labs(title = setsHeaders[i]) +
    theme(plot.title = element_text(hjust = 0.1, size=rel(0.5)), axis.title = element_text(size=rel(0.5)))
  return(p)
}


linePlotCisplatin <- function(i) {
  p <- ggplot(df) + 
    geom_line(aes(x=positions, y= df[sets[i]]/(df$GG + df$CC)), colour= "brown", size=0.1) +
    #geom_line(aes(x=positions, y=df[sets[i*4 + 3]]/df[sets[i*4 + 4]]), colour= "green", size=0.2) +
    #ylim(350, 600) +
    ylab("Damage Per 1M Mapped Reads") +
    xlab("Position Relative to TBS Center") +
    labs(title = setsHeaders[i]) +
    theme(plot.title = element_text(hjust = 0.1, size=rel(0.5)), axis.title = element_text(size=rel(0.5)))
  return(p)
}
pdf('~/SancarLab/DamageSeq/hiSeq/TFBS_individual.pdf')
multiplot(linePlotUV(1),linePlotUV(2),linePlotUV(3),linePlotUV(4),linePlotUV(5),linePlotUV(6),linePlotUV(7),linePlotUV(8), cols=4)
multiplot(linePlotUV(9),linePlotUV(10),linePlotUV(11),linePlotUV(12),linePlotUV(13),linePlotUV(14),linePlotUV(15),linePlotUV(16), cols=4)
multiplot(linePlotCisplatin(17),linePlotCisplatin(18),linePlotCisplatin(19),linePlotCisplatin(20),linePlotCisplatin(21),linePlotCisplatin(22),linePlotCisplatin(23),linePlotCisplatin(24), cols=4)
dev.off()

