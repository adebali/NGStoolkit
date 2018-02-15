library(ggplot2)
library(scales)
library(dplyr)

setwd("/proj/sancarlb/users/ogun/sancarlabutils/projects/yeast/XR/dataDir")
d <- read.table("merged_hundred.txt", header=T, row.names=NULL)
d$length <- d$end - d$start
d$nCount <- 100*d$count/d$length
d$time <- factor(d$time, levels=c("5m", "20m", "1h"))
d$category <- factor(d$category, levels=c("opposite", "same"), labels = c("TS", "NTS"))
d$product <- factor(d$product, levels=c("64PP", "CPD"), labels = c("(6-4)PP", "CPD"))
d$position <- d$position + 1

p <- ggplot(filter(d, product=="(6-4)PP" & replicate=="A")) +
  geom_tile(aes(position, y=geneNo, fill = log(count,2))) +
  scale_fill_gradient2(low = "black",
                       high = muted("blue"),
                       name="R",
                       na.value = 'white') +
  facet_grid(~time~category) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  xlab("Gene Body") +
  ylab("")

ggsave("hundred_64.pdf", p, width=6, height=8)  

p <- ggplot(filter(d, product=="CPD" & replicate=="A")) +
  geom_tile(aes(position, y=geneNo, fill = log(count,2))) +
  scale_fill_gradient2(low = "black",
                       high = muted("blue"),
                       name="R",
                       na.value = 'white') +
  facet_grid(~time~category) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  xlab("Gene Body") +
  ylab("")

ggsave("hundred_CPD.pdf", p, width=6, height=10)  