library("ggplot2")
library("reshape2")


df <- read.table(paste("/Users/ogunadebali/SancarLab/DamageSeq/hiSeq/1_8.chmm.txt", sep = ""), header = F)
headers <- names(read.table(paste("/Users/ogunadebali/SancarLab/DamageSeq/hiSeq/1_8.chmm_headers.txt", sep = ""), header = T))
names(df) <- headers
df$length <- df$end - df$start
per <- data.frame()

SUM <- sum(as.numeric(df$length))

f <- function(state) {
  return(sum(as.numeric(df[which(df["state"] == state),]$length))/SUM)
}

stateOrder = c(
  "1_Active_Promoter",
  "2_Weak_Promoter",
  "3_Poised_Promoter",
  "4_Strong_Enhancer",
  "5_Strong_Enhancer",
  "6_Weak_Enhancer",
  "7_Weak_Enhancer",
  "8_Insulator",
  "9_Txn_Transition",
  "10_Txn_Elongation",
  "11_Weak_Txn",
  "12_Repressed",
  "13_Heterochrom/lo",
  "14_Repetitive/CNV",
  "15_Repetitive/CNV")

d <- as.data.frame(matrix())

for (s in stateOrder){
  d[s] <- f(s)
}
d$V1 <- NULL
stateColors = c("#ff0000", "#ff6969", "#cf0bc6", "#faca00", "blue", "#fffc04", "black", "#0abefe", "#00b050", "red", "#99ff66", "#7f7f7f", "#c6c6c6", "red", "#c6c6c6")

e <- melt(d)
e$no = 1
ggplot(e) + geom_bar(stat="identity", aes(x=no, y=value, fill=variable)) +
  scale_fill_manual(values = stateColors)
