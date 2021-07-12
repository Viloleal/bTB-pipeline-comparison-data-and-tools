library(dplyr)
library(ggplot2)
library(remotes)
library("egg")
library(stats)
library(tidyr)
library(purrr)
library(ggpattern)
library(ggpubr)
library(grid)

#A = unfiltered, B=10bp, C=Default, , D = PE/PPE, E=PE/PPE mobile F=PPE mobile others

homop <- read.csv(file = '~/R/homoplasies/homoplasy_summary.csv')
homop$Pipeline <- factor(homop$Pipeline, levels=c("Simulation","vSNP","SNiPgenie","BovTB","MTBseq"))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
homoplot <- ggplot(homop,aes(x=Pipeline,y=Homoplasies,fill=Filter, label=Percentage)) + 
  geom_bar(position="stack", stat="identity",color="black") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size = 10),
        axis.text.x = element_text(size = 8,color="black", face = "bold"),
        axis.text.y = element_text(size = 8,color="black"),
        legend.title = element_text(size = 10, color ="black", face ="bold"),
        legend.text = element_text(size = 8,color="black"))
homoplot
