library(phytools)
library(ape)
library(ggplot2)
library(egg)

assoc <- cbind(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21",
                 "22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39",
                 "40","41","42","43","44","45"),
               c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21",
                 "22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39",
                 "40","41","42","43","44","45"))
#UNFILTERED
sim.unfilt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/sim_vcfs/core_trees/core/unfiltered/RAxML_bipartitions.sim_unfilt_core')
vsnp.unfilt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/vsnp/unfilt/RAxML_bipartitions.vsnp_unfiltered')
snpgenie.unfilt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/snpgenie/raxml_trees/unfilt/RAxML_bipartitions.snpgenie_unfiltered')
bovtb.unfilt_core <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/bovtb/trees/corrected/unfilt/corrected/RAxML_bipartitions.bovtb_unfiltered_corrected_core')
mtbseq.unfilt <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/mtbseq/raxml_trees/unfilt2/RAxML_bipartitions.mtbseq_unfiltered2')

trsim_unfilt <- drop.tip(sim.unfilt,tip = c("182","161"))
trvsnp_unfilt <- drop.tip(vsnp.unfilt,tip = c("182","161"))
trsnpgenie_unfilt <- drop.tip(snpgenie.unfilt,tip = c("182","161"))
trmtbseq_unfilt <- drop.tip(mtbseq.unfilt,tip = c("182","161"))
trbovtb.unfilt_core <- drop.tip(bovtb.unfilt_core,tip = c("182","161"))


unfilt1 <- cophylo(trsim_unfilt,trvsnp_unfilt,assoc,rotate=TRUE)
plot(unfilt1, link.lwd=2,link.col="red",fsize=0.6, scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")

unfilt2 <- cophylo(trsim_unfilt,trsnpgenie_unfilt,assoc,rotate=TRUE)
plot(unfilt2, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")

unfilt3 <- cophylo(trsim_unfilt,trmtbseq_unfilt,assoc,rotate=TRUE)
plot(unfilt3, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")

unfilt4 <- cophylo(trsim_unfilt,trbovtb.unfilt_core,assoc,rotate=TRUE)
plot(unfilt4, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")


#FILTERED
sim.tree <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/sim_vcfs/core_trees/core/10_pemobel_other/RAxML_bipartitions.sim_10_pemobel_other_filtered')
vsnp.tree <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/vsnp/10_pemobel_other_filt/RAxML_bipartitions.vsnp_10_pemobel_other_filt')
vsnpstand.tree <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/vsnp/10_standard/RAxML_bipartitions.vsnp_10_standard')
snpgenie.tree <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/snpgenie/raxml_trees/10_pemobel_other_filt/RAxML_bipartitions.snpgenie_10_pemobel_other_filt')
snpgenie_stand.tree <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/snpgenie/raxml_trees/standard/RAxML_bipartitions.snpgenie_standard')
bovtb.tree <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/bovtb/trees/corrected/10_pemobel_other_filt/corrected/RAxML_bestTree.bovtb10_pemobel_corrected_core')
bovtb_stand.tree <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/bovtb/trees/corrected/10_standard/corrected/RAxML_bipartitions.bovtb_10_standard_corrected_core')
mtbseq.tree <- read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/mtbseq/backup/RAxML_bipartitions.mtbseq_10_pemobel_other_filt2_upper')
sim_vsnp.tree <-read.tree(file = '~/Data/Wicklow/FINAL/Pipeline_results/sim_vcfs/vsnp_standard/core_fasta/RAxML_bipartitions.vsnp_standard_sim_core')
sim_genie.tree <- read.tree(file= '~/Data/Wicklow/FINAL/Pipeline_results/sim_vcfs/snpgenie_standard/core_fasta/RAxML_bipartitions.sim_standard_snpgenie_core')
sim_bovtb.tree <- read.tree(file= '~/Data/Wicklow/FINAL/Pipeline_results/sim_vcfs/bovtb_standard/core_fasta/RAxML_bipartitions.bovtb_standard_sim_core')
#drop outgroups
trsim_filt <- drop.tip(sim.tree,tip = c("182","161"))
trvsnp_filt <- drop.tip(vsnp.tree,tip = c("182","161"))
trsnpgenie_filt <- drop.tip(snpgenie.tree,tip = c("182","161"))
trmtbseq_filt <- drop.tip(mtbseq.tree,tip = c("182","161"))
trvsnp_stand <- drop.tip(vsnpstand.tree, tip = c("182","161"))
trsnpgenie_stand <- drop.tip(snpgenie_stand.tree,tip = c("182","161"))
trbovtb_filt_core <- drop.tip(bovtb.tree,tip = c("182","161"))
trbovtb_standard_core <- drop.tip(bovtb_stand.tree,tip = c("182","161"))
sim_vsnpstandard <- drop.tip(sim_vsnp.tree, tip = c("182","161"))
sim_geniestandard <- drop.tip(sim_genie.tree, tip = c("182","161"))
sim_bovtbstandard <- drop.tip(sim_bovtb.tree, tip = c("182","161"))

filtx <- cophylo(sim_vsnpstandard,trvsnp_stand,assoc,rotate=TRUE,rotate.multi=TRUE)
plot(filtx, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")

filty <- cophylo(sim_geniestandard,trsnpgenie_stand,assoc,rotate=TRUE,rotate.multi=TRUE)
plot(filty, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")

filtz <- cophylo(sim_bovtbstandard,trbovtb_standard_core,assoc,rotate=TRUE,rotate.multi=TRUE)
plot(filtz, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")

filt1 <- cophylo(trsim_filt,trvsnp_filt,assoc,rotate=TRUE,rotate.multi=TRUE)
plot(filt1, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")

filt2 <- cophylo(trsim_filt,trsnpgenie_filt,assoc,rotate=TRUE,rotate.multi=TRUE)
plot(filt2, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")

filt3 <- cophylo(trsim_filt,trmtbseq_filt,assoc,rotate=TRUE,rotate.multi=TRUE)
plot(filt3, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")

filt6 <- cophylo(trsim_filt,trvsnp_stand, assoc, rotate=TRUE,rotate.multi=TRUE)
plot(filt6, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")

filt7 <- cophylo(trvsnp_filt,trvsnp_stand, assoc, rotate=TRUE,rotate.multi=TRUE)
plot(filt7, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")

filt8 <- cophylo(trsim_filt,trsnpgenie_stand, assoc, rotate=TRUE,rotate.multi=TRUE)
plot(filt8, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")

filt9 <- cophylo(trsnpgenie_filt,trsnpgenie_stand, assoc, rotate=TRUE,rotate.multi=TRUE)
plot(filt9, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")

filt14 <- cophylo(trsim_filt,trbovtb_filt_core, assoc, rotate=TRUE,rotate.multi=TRUE)
plot(filt14, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")

filt15 <- cophylo(trsim_filt,trbovtb_standard_core, assoc, rotate=TRUE,rotate.multi=TRUE)
plot(filt15, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
plot(filt15, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")

filt16 <- cophylo(trsim_unfilt,trbovtb_standard_core, assoc, rotate=TRUE,rotate.multi=TRUE)

filt17 <- cophylo(trsim_unfilt,trsim_filt, assoc, rotate=TRUE,rotate.multi=TRUE)
plot(filt17, link.lwd=2,link.col="red",fsize=0.6,scale.bar=c(0.01,0.01))
nodelabels.cophylo(cex=0.6, frame="none",adj=c(1.1,-0.4),col="black")
nodelabels.cophylo(which="right",cex=0.6,frame="none",adj=c(-0.1,-0.4),col="black")
