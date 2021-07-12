
library(VennDiagram)
library("nVennR")
library("qpcR")
library("grImport2")
library(rsvg)

#UNFILTERED
simulated <- read.csv(file=sprintf("~/Data/Wicklow/FINAL/venn/unfiltered/simulated_positions.tsv"),header = TRUE)
genie <- read.csv(file=sprintf("~/Data/Wicklow/FINAL/venn/unfiltered/genie_positions.tsv"),header = TRUE)
vsnp <- read.csv(file=sprintf("~/Data/Wicklow/FINAL/venn/unfiltered/vsnp_positions.tsv"),header = TRUE)
mtbseq <- read.csv(file=sprintf("~/Data/Wicklow/FINAL/venn/unfiltered/new/mtbseq_positions.tsv"),header = TRUE)
bovtb <- read.csv(file=sprintf("~/Data/Wicklow/FINAL/venn/unfiltered/bovtb_positions_rem_mixed.tsv"),header = TRUE)
obj_list <- list(simulated$simulated, genie$snpgenie,vsnp$vsnp,mtbseq$mtbseq, bovtb$bovtb)

cbbPalette <- c("#999999", "#E69F00", "#0072B2", "#009E73", "#F0E442")
venn.diagram(x=obj_list,
             category.names = c("Simulated","SNiPgenie","vSNP","MTBseq","BovTB"),
             filename = sprintf('~/Data/Wicklow/FINAL/venn/unfiltered/Unfiltered_venn_final.png'),
             output=TRUE,
             euler.d=FALSE,
             scaled=FALSE,
             imagetype = "png",
             height = 2000,
             width = 2000,
             resolution = 300,
             col="black",
             lwd = 2,
             fill=cbbPalette,
             cat.col = cbbPalette,
             cat.cex = 1.4,
             cat.default.pos="outer",
             cat.fontface = "bold",
             cat.fontfamily = "serif",
             cat.pos = c(20,-15,-135,150,30),
             margin = 0.3,
             main = 'Unfiltered',
             cex = 1.3,
             fontfamiliy="serif")
#FILTERED
simulated <- read.csv(file=sprintf("~/Data/Wicklow/FINAL/venn/10_pemobel_other_filtered/simulated_positions.tsv"),header = TRUE)
genie <- read.csv(file=sprintf("~/Data/Wicklow/FINAL/venn/10_pemobel_other_filtered/genie_positions.tsv"),header = TRUE)
vsnp <- read.csv(file=sprintf("~/Data/Wicklow/FINAL/venn/10_pemobel_other_filtered/vsnp_positions.tsv"),header = TRUE)
mtbseq <- read.csv(file=sprintf("~/Data/Wicklow/FINAL/venn/10_pemobel_other_filtered/mtbseq_positions.tsv"),header = TRUE)
bovtb <- read.csv(file=sprintf("~/Data/Wicklow/FINAL/venn/10_pemobel_other_filtered/bovtb_positions_rem_mixed.tsv"),header = TRUE)
obj_list <- list(simulated$simulated, genie$snpgenie,vsnp$vsnp,mtbseq$mtbseq, bovtb$bovtb)

ab <- get.venn.partitions(obj_list, keep.elements=TRUE)
make.truth.table()
cbbPalette <- c("#999999", "#E69F00", "#0072B2", "#009E73", "#F0E442")
venn.diagram(x=obj_list,
             category.names = c("Simulated","SNiPgenie","vSNP","MTBseq","BovTB"),
             filename = sprintf('~/Data/Wicklow/FINAL/venn/10_pemobel_other_filtered/Filtered_venn_final.png'),
             output=TRUE,
             euler.d=FALSE,
             scaled=FALSE,
             imagetype = "png",
             height = 2000,
             width = 2000,
             resolution = 300,
             col="black",
             lwd = 2,
             fill=cbbPalette,
             cat.col = cbbPalette,
             cat.cex = 1.4,
             cat.default.pos="outer",
             cat.fontface = "bold",
             cat.fontfamily = "serif",
             cat.pos = c(20,-15,-135,150,30),
             margin = 0.3,
             main = 'Filtered',
             cex = 1.3,
             fontfamiliy="serif")

#FN
genie <- read.csv(file=sprintf("~/R/false_negatives/snpgenie_FN_allpositions.tsv"),header = TRUE)
vsnp <- read.csv(file=sprintf("~/R/false_negatives/vsnp_FN_allpositions.tsv"),header = TRUE)
mtbseq <- read.csv(file=sprintf("~/R/false_negatives/mtbseq_FN_allpositions.tsv"),header = TRUE)
bovtb <- read.csv(file=sprintf("~/R/false_negatives/bovtb_FN_positions.tsv"),header = TRUE)

obj_list <- list(genie$vsnp,vsnp$vsnp,mtbseq$vsnp, bovtb$vsnp)
cbbPalette <- c("#E69F00", "#0072B2", "#009E73", "#F0E442")
venn.diagram(x=obj_list,
             category.names = c("SNiPgenie","vSNP","MTBseq","BovTB"),
             filename = sprintf('~/R/false_negatives/FN_venn_final.png'),
             output=TRUE,
             euler.d=FALSE,
             scaled=FALSE,
             imagetype = "png",
             height = 2000,
             width = 2000,
             resolution = 300,
             col="black",
             lwd = 2,
             fill=cbbPalette,
             cat.cex = 1.4,
             cat.default.pos="outer",
             cat.fontface = "bold",
             cat.fontfamily = "serif",
             margin = 0.3,
             main = 'Unfiltered',
             cex = 1.3,
             fontfamiliy="serif")
