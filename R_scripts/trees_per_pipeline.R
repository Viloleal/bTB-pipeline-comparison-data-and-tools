library(treespace)
library(phytools)
library("adegenet")
library("adegraphics")
library("rgl")
library(ggplot2)
library(rgl)
library(funData)
library("RColorBrewer")

sim_unfilt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/sim_vcfs/core_trees/core/unfiltered/RAxML_bootstrap.sim_unfilt_core')
sim_10_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/sim_vcfs/core_trees/core/10_filt/RAxML_bootstrap.sim_10_filt')
sim_10_pe_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/sim_vcfs/core_trees/core/10_pe_filt/RAxML_bootstrap.sim_10_pe_filt')
sim_10_pemobel_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/sim_vcfs/core_trees/core/10_pemobel_filt/RAxML_bootstrap.sim_10_pemobel_filt')
sim_10_pemobel_other_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/sim_vcfs/core_trees/core/10_pemobel_other/RAxML_bootstrap.sim_10_pemobel_other_filtered')

sim_unfilt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/sim_vcfs/core_trees/core/unfiltered/RAxML_bestTree.sim_unfilt_core')
sim_10_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/sim_vcfs/core_trees/core/10_filt/RAxML_bestTree.sim_10_filt')
sim_10_pe_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/sim_vcfs/core_trees/core/10_pe_filt/RAxML_bestTree.sim_10_pe_filt')
sim_10_pemobel_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/sim_vcfs/core_trees/core/10_pemobel_filt/RAxML_bestTree.sim_10_pemobel_filt')
sim_10_pemobel_other_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/sim_vcfs/core_trees/core/10_pemobel_other/RAxML_bestTree.sim_10_pemobel_other_filtered')
all_trees <- c(sim_unfilt, sim_10_filt, sim_10_pe_filt, sim_10_pemobel_filt, sim_10_pemobel_other_filt,
               sim_unfilt_b, sim_10_filt_b, sim_10_pe_filt_b, sim_10_pemobel_filt_b, sim_10_pemobel_other_filt_b)
class(all_trees) <- "multiPhylo"

names(all_trees)[1:100] <- paste0("sim_unfilt",1:100)
names(all_trees)[101:200] <- paste0("sim_10_filt",1:100)
names(all_trees)[201:300] <- paste0("sim_10_pe_filt",1:100)
names(all_trees)[301:400] <- paste0("sim_10_pemobel_filt",1:100)
names(all_trees)[401:500] <- paste0("sim_10_pemobel_other",1:100)
names(all_trees)[501] <- paste0("sim_unfilt_best")
names(all_trees)[502] <- paste0("sim_10_filt_best")
names(all_trees)[503] <- paste0("sim_10_pe_filt_best")
names(all_trees)[504] <- paste0("sim_10_pemobel_filt_best")
names(all_trees)[505] <- paste0("sim_10_pemobel_other_filt_best")

Dtype <- c(rep("Unfiltered",100),rep("10bp window",100),rep("10bp + PE/PPE family proteins",100),
           rep("10bp + PE/PPE + mobile elements",100), rep("10bp + PE/PPE + mobile elements + other",100),
           rep("Unfiltered",1),rep("10bp window",1),rep("10bp + PE/PPE family proteins",1),
           rep("10bp + PE/PPE + mobile elements",1), rep("10bp + PE/PPE + mobile elements + other",1))
           
Dcols <- c("#56B4E9", "#009E73", "#E69F00", "#0072B2", "#D55E00")
      
Dscape <- treespace(all_trees,method='RF', nf=10)
symbs <- c(rep("Bootstrap replicate",500),rep("Best tree",5)) 
plotGrovesD3(Dscape$pco,
             fixed=TRUE,
             groups=Dtype, 
             tooltip_text = names(all_trees), # add the tree names as tooltip text
             colors=Dcols,
             col_lab="Pipeline",
             legend_width=150,
             point_size = 30,
             symbol_var = symbs,
             ellipses=TRUE)

#Bovtb
bovtb_unfilt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/bovtb/trees/corrected/unfilt/corrected/RAxML_bootstrap.bovtb_unfiltered_corrected_core')
bovtb_10_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/bovtb/trees/corrected/10_filt/corrected/RAxML_bootstrap.bovtb_10_filt_corrected_core')
bovtb_10_pe_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/bovtb/trees/corrected/10_pe_filt/corrected/RAxML_bootstrap.bovtb_10_pe_filt_corrected_core')
bovtb_10_pemobel_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/bovtb/trees/corrected/10_pemobel_filt/corrected/RAxML_bootstrap.bovtb_10_pemobel_filt_corrected_core')
bovtb_10_pemobel_other_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/bovtb/trees/corrected/10_pemobel_other_filt/corrected/RAxML_bootstrap.bovtb10_pemobel_corrected_core')
bovtb_10_default_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/bovtb/trees/corrected/10_standard/corrected/RAxML_bootstrap.bovtb_10_standard_corrected_core')

bovtb_unfilt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/bovtb/trees/corrected/unfilt/corrected/RAxML_bestTree.bovtb_unfiltered_corrected_core')
bovtb_10_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/bovtb/trees/corrected/10_filt/corrected/RAxML_bestTree.bovtb_10_filt_corrected_core')
bovtb_10_pe_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/bovtb/trees/corrected/10_pe_filt/corrected/RAxML_bestTree.bovtb_10_pe_filt_corrected_core')
bovtb_10_pemobel_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/bovtb/trees/corrected/10_pemobel_filt/corrected/RAxML_bestTree.bovtb_10_pemobel_filt_corrected_core')
bovtb_10_pemobel_other_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/bovtb/trees/corrected/10_pemobel_other_filt/corrected/RAxML_bestTree.bovtb10_pemobel_corrected_core')
bovtb_10_default_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/bovtb/trees/corrected/10_standard/corrected/RAxML_bestTree.bovtb_10_standard_corrected_core')
all_trees <- c(sim_unfilt, sim_10_filt, sim_10_pe_filt, sim_10_pemobel_filt, sim_10_pemobel_other_filt,
               sim_unfilt_b, sim_10_filt_b, sim_10_pe_filt_b, sim_10_pemobel_filt_b, sim_10_pemobel_other_filt_b,
               bovtb_unfilt, bovtb_10_filt, bovtb_10_pe_filt, bovtb_10_pemobel_filt, bovtb_10_pemobel_other_filt,bovtb_10_default_filt,
               bovtb_unfilt_b, bovtb_10_filt_b, bovtb_10_pe_filt_b, bovtb_10_pemobel_filt_b, bovtb_10_pemobel_other_filt_b,bovtb_10_default_filt_b)
class(all_trees) <- "multiPhylo"
names(all_trees)[1:100] <- paste0("sim_unfilt",1:100)
names(all_trees)[101:200] <- paste0("sim_10_filt",1:100)
names(all_trees)[201:300] <- paste0("sim_10_pe_filt",1:100)
names(all_trees)[301:400] <- paste0("sim_10_pemobel_filt",1:100)
names(all_trees)[401:500] <- paste0("sim_10_pemobel_other",1:100)
names(all_trees)[501] <- paste0("sim_unfilt_best")
names(all_trees)[502] <- paste0("sim_10_filt_best")
names(all_trees)[503] <- paste0("sim_10_pe_filt_best")
names(all_trees)[504] <- paste0("sim_10_pemobel_filt_best")
names(all_trees)[505] <- paste0("sim_10_pemobel_other_filt_best")
names(all_trees)[506:605] <- paste0("bovtb_unfilt",1:100)
names(all_trees)[606:705] <- paste0("bovtb_10_filt",1:100)
names(all_trees)[706:805] <- paste0("bovtb_10_pe_filt",1:100)
names(all_trees)[806:905] <- paste0("bovtb_10_pemobel_filt",1:100)
names(all_trees)[906:1005] <- paste0("bovtb_10_pemobel_other",1:100)
names(all_trees)[1006:1105] <- paste0("bovtb_10_default",1:100)
names(all_trees)[1106] <- paste0("bovtb_unfilt_best")
names(all_trees)[1107] <- paste0("bovtb_10_filt_best")
names(all_trees)[1108] <- paste0("bovtb_10_pe_filt_best")
names(all_trees)[1109] <- paste0("bovtb_10_pemobel_filt_best")
names(all_trees)[1110] <- paste0("bovtb_10_pemobel_other_filt_best")
names(all_trees)[1111] <- paste0("bovtb_10_default_best")
Dtype <- c(rep("Unfiltered",100),rep("10bp window",100),rep("10bp + PE/PPE family proteins",100),
           rep("10bp + PE/PPE + mobile elements",100), rep("10bp + PE/PPE + mobile elements + other",100),
           rep("Best tree",5),
           rep("Unfiltered",100),rep("10bp window",100),rep("10bp + PE/PPE family proteins",100),
           rep("10bp + PE/PPE + mobile elements",100), rep("10bp + PE/PPE + mobile elements + other",100),
           rep("10bp + default filters",100),
           rep("Best tree",6))
Dcols <- c("#999999","#56B4E9", "#009E73", "#E69F00", "FF624E", "#D55E00", "#CC79A7")

symbs <- c(rep("Simulated",505),
           rep("BovTB",606))
Dscape <- treespace(all_trees,method='RF', nf=12)
plotGrovesD3(Dscape$pco,
             fixed=TRUE,
             groups=Dtype, 
             tooltip_text = names(all_trees), # add the tree names as tooltip text
             colors=Dcols,
             col_lab="Pipeline",
             legend_width=150,
             point_size = 30,
             symbol_var = symbs,
             symbols = c("Simulated"="circle","BovTB"="triangle"),
             ellipses=TRUE)
#MTBseq
mtbseq_unfilt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/mtbseq/raxml_trees/unfilt2/RAxML_bootstrap.mtbseq_unfiltered2')
mtbseq_10_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/mtbseq/raxml_trees/10_filt//RAxML_bootstrap.mtbseq_10_filt2')
mtbseq_10_pe_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/mtbseq/raxml_trees/10_pe_filt/RAxML_bootstrap.mtbseq_10_pe_filt2')
mtbseq_10_pemobel_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/mtbseq/raxml_trees/10_pemobel_filt/RAxML_bootstrap.mtbseq_10_pemobel_filt2')
mtbseq_10_pemobel_other_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/mtbseq/raxml_trees/10_pemobel_other_filt2/RAxML_bootstrap.mtbseq_10_pemobel_other_filt2')

mtbseq_unfilt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/mtbseq/raxml_trees/unfilt2/RAxML_bestTree.mtbseq_unfiltered2')
mtbseq_10_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/mtbseq/raxml_trees/10_filt//RAxML_bestTree.mtbseq_10_filt2')
mtbseq_10_pe_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/mtbseq/raxml_trees/10_pe_filt/RAxML_bestTree.mtbseq_10_pe_filt2')
mtbseq_10_pemobel_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/mtbseq/raxml_trees/10_pemobel_filt/RAxML_bestTree.mtbseq_10_pemobel_filt2')
mtbseq_10_pemobel_other_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/mtbseq/raxml_trees/10_pemobel_other_filt2/RAxML_bestTree.mtbseq_10_pemobel_other_filt2')

all_trees <- all_trees <- c(sim_unfilt, sim_10_filt, sim_10_pe_filt, sim_10_pemobel_filt, sim_10_pemobel_other_filt,
                            sim_unfilt_b, sim_10_filt_b, sim_10_pe_filt_b, sim_10_pemobel_filt_b, sim_10_pemobel_other_filt_b,
                            mtbseq_unfilt, mtbseq_10_filt, mtbseq_10_pe_filt, mtbseq_10_pemobel_filt, mtbseq_10_pemobel_other_filt,
                            mtbseq_unfilt_b, mtbseq_10_filt_b, mtbseq_10_pe_filt_b, mtbseq_10_pemobel_filt_b, mtbseq_10_pemobel_other_filt_b)
class(all_trees) <- "multiPhylo"
names(all_trees)[1:100] <- paste0("sim_unfilt",1:100)
names(all_trees)[101:200] <- paste0("sim_10_filt",1:100)
names(all_trees)[201:300] <- paste0("sim_10_pe_filt",1:100)
names(all_trees)[301:400] <- paste0("sim_10_pemobel_filt",1:100)
names(all_trees)[401:500] <- paste0("sim_10_pemobel_other",1:100)
names(all_trees)[501] <- paste0("sim_unfilt_best")
names(all_trees)[502] <- paste0("sim_10_filt_best")
names(all_trees)[503] <- paste0("sim_10_pe_filt_best")
names(all_trees)[504] <- paste0("sim_10_pemobel_filt_best")
names(all_trees)[505] <- paste0("sim_10_pemobel_other_filt_best")
names(all_trees)[506:605] <- paste0("mtbseq_unfilt",1:100)
names(all_trees)[606:705] <- paste0("mtbseq_10_filt",1:100)
names(all_trees)[706:805] <- paste0("mtbseq_10_pe_filt",1:100)
names(all_trees)[806:905] <- paste0("mtbseq_10_pemobel_filt",1:100)
names(all_trees)[906:1005] <- paste0("mtbseq_10_pemobel_other",1:100)
names(all_trees)[1006] <- paste0("mtbseq_unfilt_best")
names(all_trees)[1007] <- paste0("mtbseq_10_filt_best")
names(all_trees)[1008] <- paste0("mtbseq_10_pe_filt_best")
names(all_trees)[1009] <- paste0("mtbseq_10_pemobel_filt_best")
names(all_trees)[1010] <- paste0("mtbseq_10_pemobel_other_filt_best")
Dtype <- c(rep("Unfiltered",100),rep("10bp window",100),rep("10bp + PE/PPE family proteins",100),
           rep("10bp + PE/PPE + mobile elements",100), rep("10bp + PE/PPE + mobile elements + other",100),
           rep("Best tree",5),
           rep("Unfiltered",100),rep("10bp window",100),rep("10bp + PE/PPE family proteins",100),
           rep("10bp + PE/PPE + mobile elements",100), rep("10bp + PE/PPE + mobile elements + other",100),
           rep("Best tree",5))
Dcols <- c("#999999","#56B4E9", "#009E73", "#E69F00", "FF624E", "#D55E00", "#CC79A7")

symbs <- c(rep("Simulated",505),
           rep("MTBseq",505))
Dscape <- treespace(all_trees,method='RF', nf=10)

plotGrovesD3(Dscape$pco,
             fixed=TRUE,
             groups=Dtype, 
             tooltip_text = names(all_trees), # add the tree names as tooltip text
             colors=Dcols,
             col_lab="Pipeline",
             legend_width=150,
             point_size = 30,
             symbol_var = symbs,
             symbols = c("Simulated"="circle","MTBseq"="triangle"),
             ellipses=TRUE)
#vSNP
vsnp_unfilt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/vsnp/unfilt/RAxML_bootstrap.vsnp_unfiltered')
vsnp_10_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/vsnp/10_filt/RAxML_bootstrap.vsnp_10_filt')
vsnp_10_pe_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/vsnp/10_pe_filt/RAxML_bootstrap.vsnp_10_pe_filt')
vsnp_10_pemobel_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/vsnp/10_pemobel_filt/RAxML_bootstrap.vsnp_10_pemobel_filt')
vsnp_10_pemobel_other_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/vsnp/10_pemobel_other_filt/RAxML_bootstrap.vsnp_10_pemobel_other_filt')
vsnp_10_default_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/vsnp/10_standard/RAxML_bootstrap.vsnp_10_standard')

vsnp_unfilt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/vsnp/unfilt/RAxML_bestTree.vsnp_unfiltered')
vsnp_10_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/vsnp/10_filt/RAxML_bestTree.vsnp_10_filt')
vsnp_10_pe_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/vsnp/10_pe_filt/RAxML_bestTree.vsnp_10_pe_filt')
vsnp_10_pemobel_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/vsnp/10_pemobel_filt/RAxML_bestTree.vsnp_10_pemobel_filt')
vsnp_10_pemobel_other_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/vsnp/10_pemobel_other_filt/RAxML_bestTree.vsnp_10_pemobel_other_filt')
vsnp_10_default_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/vsnp/10_standard/RAxML_bestTree.vsnp_10_standard')

all_trees <- c(sim_unfilt, sim_10_filt, sim_10_pe_filt, sim_10_pemobel_filt, sim_10_pemobel_other_filt,
               sim_unfilt_b, sim_10_filt_b, sim_10_pe_filt_b, sim_10_pemobel_filt_b, sim_10_pemobel_other_filt_b,
               vsnp_unfilt, vsnp_10_filt, vsnp_10_pe_filt, vsnp_10_pemobel_filt, vsnp_10_pemobel_other_filt,vsnp_10_default_filt,
               vsnp_unfilt_b, vsnp_10_filt_b, vsnp_10_pe_filt_b, vsnp_10_pemobel_filt_b, vsnp_10_pemobel_other_filt_b,vsnp_10_default_filt_b)
class(all_trees) <- "multiPhylo"
names(all_trees)[1:100] <- paste0("sim_unfilt",1:100)
names(all_trees)[101:200] <- paste0("sim_10_filt",1:100)
names(all_trees)[201:300] <- paste0("sim_10_pe_filt",1:100)
names(all_trees)[301:400] <- paste0("sim_10_pemobel_filt",1:100)
names(all_trees)[401:500] <- paste0("sim_10_pemobel_other",1:100)
names(all_trees)[501] <- paste0("sim_unfilt_best")
names(all_trees)[502] <- paste0("sim_10_filt_best")
names(all_trees)[503] <- paste0("sim_10_pe_filt_best")
names(all_trees)[504] <- paste0("sim_10_pemobel_filt_best")
names(all_trees)[505] <- paste0("sim_10_pemobel_other_filt_best")
names(all_trees)[506:605] <- paste0("vsnp_unfilt",1:100)
names(all_trees)[606:705] <- paste0("vsnp_10_filt",1:100)
names(all_trees)[706:805] <- paste0("vsnp_10_pe_filt",1:100)
names(all_trees)[806:905] <- paste0("vsnp_10_pemobel_filt",1:100)
names(all_trees)[906:1005] <- paste0("vsnp_10_pemobel_other",1:100)
names(all_trees)[1006:1105] <- paste0("vsnp_10_default",1:100)
names(all_trees)[1106] <- paste0("vsnp_unfilt_best")
names(all_trees)[1107] <- paste0("vsnp_10_filt_best")
names(all_trees)[1108] <- paste0("vsnp_10_pe_filt_best")
names(all_trees)[1109] <- paste0("vsnp_10_pemobel_filt_best")
names(all_trees)[1110] <- paste0("vsnp_10_pemobel_other_filt_best")
names(all_trees)[1111] <- paste0("vsnp_10_default_best")
Dtype <- c(rep("Unfiltered",100),rep("10bp window",100),rep("10bp + PE/PPE family proteins",100),
           rep("10bp + PE/PPE + mobile elements",100), rep("10bp + PE/PPE + mobile elements + other",100),
           rep("Best tree",5),
           rep("Unfiltered",100),rep("10bp window",100),rep("10bp + PE/PPE family proteins",100),
           rep("10bp + PE/PPE + mobile elements",100), rep("10bp + PE/PPE + mobile elements + other",100),
           rep("10bp + default filters",100),
           rep("Best tree",6))
Dcols <- c("#999999","#56B4E9", "#009E73", "#E69F00", "FF624E", "#D55E00", "#CC79A7")

symbs <- c(rep("Simulated",505),
           rep("vSNP",606))
Dscape <- treespace(all_trees,method='RF', nf=12)
plotGrovesD3(Dscape$pco,
             fixed=TRUE,
             groups=Dtype, 
             tooltip_text = names(all_trees), # add the tree names as tooltip text
             colors=Dcols,
             col_lab="Pipeline",
             legend_width=150,
             point_size = 30,
             symbol_var = symbs,
             symbols = c("Simulated"="circle","vSNP"="triangle"),
             ellipses=TRUE)
#SNiPgenie
snpgenie_unfilt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/snpgenie/raxml_trees/unfilt/RAxML_bootstrap.snpgenie_unfiltered')
snpgenie_10_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/snpgenie/raxml_trees/10_filt//RAxML_bootstrap.snpgenie_10_filt')
snpgenie_10_pe_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/snpgenie/raxml_trees/10_pe_filt/RAxML_bootstrap.snpgenie_10_pe_filt')
snpgenie_10_pemobel_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/snpgenie/raxml_trees/10_pemobel_filt/RAxML_bootstrap.snpgenie_10_pemobel_filt')
snpgenie_10_pemobel_other_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/snpgenie/raxml_trees/10_pemobel_other_filt/RAxML_bootstrap.snpgenie_10_pemobel_other_filt')
snpgenie_10_default_filt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/snpgenie/raxml_trees/standard/RAxML_bootstrap.snpgenie_standard')

snpgenie_unfilt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/snpgenie/raxml_trees/unfilt/RAxML_bestTree.snpgenie_unfiltered')
snpgenie_10_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/snpgenie/raxml_trees/10_filt//RAxML_bestTree.snpgenie_10_filt')
snpgenie_10_pe_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/snpgenie/raxml_trees/10_pe_filt/RAxML_bestTree.snpgenie_10_pe_filt')
snpgenie_10_pemobel_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/snpgenie/raxml_trees/10_pemobel_filt/RAxML_bestTree.snpgenie_10_pemobel_filt')
snpgenie_10_pemobel_other_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/snpgenie/raxml_trees/10_pemobel_other_filt/RAxML_bestTree.snpgenie_10_pemobel_other_filt')
snpgenie_10_default_filt_b <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/snpgenie/raxml_trees/standard/RAxML_bestTree.snpgenie_standard')

all_trees <- c(sim_unfilt, sim_10_filt, sim_10_pe_filt, sim_10_pemobel_filt, sim_10_pemobel_other_filt,
               sim_unfilt_b, sim_10_filt_b, sim_10_pe_filt_b, sim_10_pemobel_filt_b, sim_10_pemobel_other_filt_b,
               snpgenie_unfilt, snpgenie_10_filt, snpgenie_10_pe_filt, snpgenie_10_pemobel_filt, snpgenie_10_pemobel_other_filt,snpgenie_10_default_filt,
               snpgenie_unfilt_b, snpgenie_10_filt_b, snpgenie_10_pe_filt_b, snpgenie_10_pemobel_filt_b, snpgenie_10_pemobel_other_filt_b,snpgenie_10_default_filt_b)
class(all_trees) <- "multiPhylo"
names(all_trees)[1:100] <- paste0("sim_unfilt",1:100)
names(all_trees)[101:200] <- paste0("sim_10_filt",1:100)
names(all_trees)[201:300] <- paste0("sim_10_pe_filt",1:100)
names(all_trees)[301:400] <- paste0("sim_10_pemobel_filt",1:100)
names(all_trees)[401:500] <- paste0("sim_10_pemobel_other",1:100)
names(all_trees)[501] <- paste0("sim_unfilt_best")
names(all_trees)[502] <- paste0("sim_10_filt_best")
names(all_trees)[503] <- paste0("sim_10_pe_filt_best")
names(all_trees)[504] <- paste0("sim_10_pemobel_filt_best")
names(all_trees)[505] <- paste0("sim_10_pemobel_other_filt_best")
names(all_trees)[506:605] <- paste0("snpgenie_unfilt",1:100)
names(all_trees)[606:705] <- paste0("snpgenie_10_filt",1:100)
names(all_trees)[706:805] <- paste0("snpgenie_10_pe_filt",1:100)
names(all_trees)[806:905] <- paste0("snpgenie_10_pemobel_filt",1:100)
names(all_trees)[906:1005] <- paste0("snpgenie_10_pemobel_other",1:100)
names(all_trees)[1006:1105] <- paste0("snpgenie_10_default",1:100)
names(all_trees)[1106] <- paste0("snpgenie_unfilt_best")
names(all_trees)[1107] <- paste0("snpgenie_10_filt_best")
names(all_trees)[1108] <- paste0("snpgenie_10_pe_filt_best")
names(all_trees)[1109] <- paste0("snpgenie_10_pemobel_filt_best")
names(all_trees)[1110] <- paste0("snpgenie_10_pemobel_other_filt_best")
names(all_trees)[1111] <- paste0("snpgenie_10_default_best")
Dtype <- c(rep("Unfiltered",100),rep("10bp window",100),rep("10bp + PE/PPE family proteins",100),
           rep("10bp + PE/PPE + mobile elements",100), rep("10bp + PE/PPE + mobile elements + other",100),
           rep("Best tree",5),
           rep("Unfiltered",100),rep("10bp window",100),rep("10bp + PE/PPE family proteins",100),
           rep("10bp + PE/PPE + mobile elements",100), rep("10bp + PE/PPE + mobile elements + other",100),
           rep("10bp + default filters",100),
           rep("Best tree",6))
Dcols <- c("#999999","#56B4E9", "#009E73", "#E69F00", "FF624E", "#D55E00", "#CC79A7")

symbs <- c(rep("Simulated",505),
           rep("snpgenie",606))
Dscape <- treespace(all_trees,method='RF', nf=13)
plotGrovesD3(Dscape$pco,
             fixed=TRUE,
             groups=Dtype, 
             tooltip_text = names(all_trees), # add the tree names as tooltip text
             colors=Dcols,
             col_lab="Pipeline",
             legend_width=150,
             point_size = 30,
             symbol_var = symbs,
             symbols = c("Simulated"="circle","snpgenie"="triangle"),
             ellipses=TRUE)

