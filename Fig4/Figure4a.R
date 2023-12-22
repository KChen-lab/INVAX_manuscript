###################################Figure 4a and 4e - GES score calculated by GSDensity ########################

library(Seurat)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(wesanderson)
library(rstatix)
#################################################Figure 4a ########################################


#GSDensity score can be calculated using GSDensity_INVAX.R
gsd=read.csv('./GSDensity_score.csv')#Each row represents a cell, in the same order for input CD8 seurat object

#check if colnames of gsd matches
identical(colnames(gsd),c('checkpoint','cytotoxic','tissue_memory'))


#This code needs meta data for CD8 T cells including cell barcode, response, patient id, timepoint; raw data can be downloaded from EGAC50000000088
cd8=readRDS('./CD8.rds')
meta=cd8@meta.data

meta=cbind(meta,gsd)
scores=colnames(gsd)
colnames(meta)[6]='patient'
d=data.frame(patient=unique(meta$patient ))
for(i in 1:nrow(d)){
  d[i,'Response']=meta$Response[which(meta$patient==unique(meta$patient)[i])][1]
  d[i,paste0('pre_',scores)]=colMeans(meta[which(meta$Timepoint==1&meta$patient==unique(meta$patient)[i]),scores])
  d[i,paste0('post_',scores)]=colMeans(meta[which(meta$Timepoint==2&meta$patient==unique(meta$patient)[i]),scores])
}
d=d[which(!is.na(d$Response)),]
d=d[which(d$patient!='INVAX033'),]

d$Response=factor(d$Response,levels=c('MHR','nonMHR'))
p1=ggboxplot(d,'Response','pre_cytotoxic',add = 'point',fill='Response',palette = c("darkblue",  "skyblue"))+stat_compare_means()+
  xlab('')+ggtitle('Cytotoxicity Score')+ylab('Gene set score')+theme(legend.position = 'none')
p2=ggboxplot(d,'Response','pre_checkpoint',add = 'point',fill='Response',palette = c("darkblue",  "skyblue"))+stat_compare_means()+
xlab('')+ggtitle('Immune checkpoint Score')+ylab('Gene set score')+theme(legend.position = 'none')

p3=ggboxplot(d,'Response','pre_tissue_memory',add = 'point',fill='Response',palette = c("darkblue",  "skyblue"))+stat_compare_means()+
  xlab('')+ggtitle('Tissue resident memory Score')+ylab('Gene set score')+theme(legend.position = 'none')
p1+p2+p3

stat.trm=wilcox_test(pre_tissue_memory~Response,data=d)
stat.check=wilcox_test(pre_checkpoint~Response,data=d)
stat.cyto=wilcox_test(pre_cytotoxic~Response,data=d)

#################################Figure 4e##################################
pal <- wes_palette("Zissou1", 100, type = "continuous")

df=data.frame(cluster=levels(meta$CellSubset))
scores=colnames(gsd)
for(i in 1:nrow(df)){
  df[i,'Abundance']=length(which(meta$CellSubset==df$cluster[i]))/nrow(meta)
  df[i,scores]= colMeans(meta[which(meta$CellSubset==df$cluster[i]),scores])
}

ggscatter(df,"cytotoxic","tissue_memory" ,color='checkpoint',
          size='Abundance')+geom_text_repel(aes(label=cluster))+scale_color_gradientn(colours = pal,name='Checkpoint\n score')+
  theme(legend.position="right")+xlab('Cytotoxicity Score')+ylab('Tissue resident Score')

df$Abundance=df$Abundance*100





