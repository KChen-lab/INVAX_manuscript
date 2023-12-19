####################################################Figure3f and 3g###########################################

library(Seurat)
library(ggplot2)
library(AUCell)
library(GSEABase)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
nk=readRDS('nk.rds')#This code needs NK seurat object

meta_nk=nk@meta.data
nk_cd16=subset(nk,subset=FCGR3A>0)
meta_nk$FCGR3A_exp='negative'
meta_nk[colnames(nk_cd16),'FCGR3A_exp']='positive'
meta_nk$FCGR3A_exp=factor(meta_nk$FCGR3A_exp,levels = c('positive','negative'))
nk@meta.data=meta_nk
#Figure 3f
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

p1=VlnPlot(nk,features = c('Cytotoxicity'),group.by = 'FCGR3A_exp',pt.size = 0,
           cols = rep('gold',3))+stat_compare_means()&
  stat_summary(fun.data=data_summary)&NoLegend()&xlab('')
p2=VlnPlot(nk,features = c('Stimulatory'),group.by = 'FCGR3A_exp',pt.size = 0,
           cols = rep('gold',3))+stat_compare_means()&
  stat_summary(fun.data=data_summary)&NoLegend()&xlab('')
p3=VlnPlot(nk,features = c('Inhibitory'),group.by = 'FCGR3A_exp',pt.size = 0,
           cols = rep('gold',3))+stat_compare_means()&
  stat_summary(fun.data=data_summary)&NoLegend()&xlab('')
p1|p2|p3
ggsave('~/Documents/ChenLab/Projects/HeadNeck/final.figures/Figure 3/NK_CD16_violin.pdf',width = 7,height = 4)

cyto.1=data_summary(nk$Cytotoxicity[which(nk$FCGR3A_exp=='positive')])
cyto.2=data_summary(nk$Cytotoxicity[which(nk$FCGR3A_exp=='negative')])
cyto=wilcox_test(Cytotoxicity~FCGR3A_exp,data=meta_nk)
stim.1=data_summary(nk$Stimulatory[which(nk$FCGR3A_exp=='positive')])
stim.2=data_summary(nk$Stimulatory[which(nk$FCGR3A_exp=='negative')])
stim=wilcox_test(Stimulatory~FCGR3A_exp,data=meta_nk)
inhi.1=data_summary(nk$Inhibitory[which(nk$FCGR3A_exp=='positive')])
inhi.2=data_summary(nk$Inhibitory[which(nk$FCGR3A_exp=='negative')])
inhi=wilcox_test(Inhibitory~FCGR3A_exp,data=meta_nk)

stats=rbind(cyto.1,cyto.2,stim.1,stim.2,inhi.1,inhi.2)
rownames(stats)=c('Cytotoxcity.Positive','Cytotoxcity.Negative',
                  'Stimulatory.Positive','Stimulatory.Negative',
                  'Inhibitory.Positive','Inhibitory.Negative')
stats.1=rbind(cyto,stim,inhi)

#Figure 3g
FeaturePlot(nk,features = c('Stimulatory','Inhibitory','Cytotoxicity'),ncol = 3,cols=rev(brewer.pal(11,'Spectral')))&NoAxes()

#######################Optional; if want to recalculate the AUCell score############################
#####AUCell score
exprMatrix <- GetAssayData(nk,slot = 'count')
cells_rankings<-AUCell_buildRankings(exprMatrix)

#Cyrotoxic score
geneset=c('PRF1','GZMH','GZMA','GZMB','GZMM','FASLG','GNLY','NKG7')
geneSets <- GeneSet(geneset, setName='Cytotoxicity')
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
nk$Cytotoxicity=as.numeric(getAUC(cells_AUC)['Cytotoxicity',]) 

#Inhibitory score
geneset=c('KIR2DL3','KIR2DL4','KLRC1','TIGIT','CD96')
geneSets <- GeneSet(geneset, setName='Inhibitory receptor')
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
nk$Inhibitory=as.numeric(getAUC(cells_AUC)['Inhibitory receptor',]) 

#Stimulatoryscore
geneset=c('FCGR3A','NCR3','KLRF1','KLRK1','KLRC2')
geneSets <- GeneSet(geneset, setName='Stimulatory receptor')
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
nk$Stimulatory=as.numeric(getAUC(cells_AUC)['Stimulatory receptor',]) 

