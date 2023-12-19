###############################################Figure 3d######################################

library(Seurat)
library(ggplot2)
library(rstatix)
library(ggpubr)
nk=readRDS('./nk.rds') #This code needs NK seurat object

meta_nk=nk@meta.data
nk_cd16=subset(nk,subset=FCGR3A>0)
meta_nk$FCGR3A_exp='None'
meta_nk[colnames(nk_cd16),'FCGR3A_exp']='CD16'
nk@meta.data=meta_nk


cell_type_proportion <- function(metadata,paired){
  stat=table(metadata$Patient,metadata$CellSubset,metadata$Timepoint)
  pre=as.matrix(stat[,,1])
  pre=pre/rowSums(pre)
  post=as.matrix(stat[,,2])
  post=post/rowSums(post)
  pre1=as.data.frame(pre)
  post1=as.data.frame(post)
  pre1$Timepoint='Baseline'
  post1$Timepoint='OnTreatment'
  df=rbind(pre1,post1)
  colnames(df)=c('Patient','CellType','Proportion','Timepoint')
  df$Proportion=df$Proportion*100
  tt=table(metadata$Patient,metadata$Response)
  tt=as.data.frame(tt)
  tt=tt[which(tt$Freq!=0),]
  rownames(tt)=tt$Var1
  df$Response=tt[as.character(df$Patient),'Var2']
  df=df[which(df$Patient!='INVAX026'),]
  return(df)
}

meta_nk$CD16=0
meta_nk[WhichCells(subset(nk,subset=FCGR3A>0)),'CD16']='CD16 positive'
meta_nk$CellSubset=meta_nk$CD16
df=cell_type_proportion(meta_nk,0)
colnames(df)
df_plot=df[which(df$Timepoint=='Baseline'&df$CellType=='CD16 positive'&!is.na(df$Response)&!is.na(df$Proportion)),]
ggplot(df_plot,aes(x=Response,y=Proportion,fill=Response))+geom_boxplot()+geom_point(size=2)+
  stat_compare_means()+scale_fill_manual(values=c('royalblue','lightblue'))+theme_classic(base_size = 15)+ylab('% Total baseline NK cells')+xlab('')+
  theme(legend.position = 'none',axis.text = element_text(color='black'))+ggtitle('Baseline FCGR3A+ NK cells')
wilcox_test(Proportion~Response,data=df_plot,conf.level  = 0.95)



