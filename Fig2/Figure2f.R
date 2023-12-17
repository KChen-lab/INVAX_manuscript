#############################################Figure2f boxplots comparing CD8/Treg ratio between response group#############

library(readxl)
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)


stat1=read_xlsx('./Supplementary Tables 6-13_v3.xlsx',
                sheet='Supplementary Table 6',range = 'A3:N73')
response=read_xlsx('../DataShare/response.xlsx')
response=as.data.frame(response)

stat1[,4:14]=stat1[,4:14]/stat1$`Total No. Cells after QC`
stat1=stat1[,-3]
colnames(stat1)
stat1=melt(stat1,id.vars = c('Patient ID','Timepoint'),variable.name ='CellType')
colnames(stat1)=c('Patient','Timepoint','CellType','Proportion')
stat1$Proportion=stat1$Proportion*100
stat1$Response=as.character(response[match(stat1$Patient,response$Patient),'Response'])
stat1$Response=factor(stat1$Response,levels=c('MHR','nMHR'))
stat1$CellType=factor(stat1$CellType,levels = c('Malignant','Epithelial','Endothelial','Fibroblast','NK','Myeloid','B','Plasma','CD8','CD4','Treg')) 

stat1$Timepoint[which(stat1$Timepoint=='OnTreatment')]='OnInduction'


df=cell_type_proportion(meta,3)
df_ratio=data.frame(
  CD8toTreg=stat1$Proportion[which(stat1$CellType=='CD8')]/stat1$Proportion[which(stat1$CellType=='Treg')],
  Timepoint=stat1$Timepoint[which(stat1$CellType=='CD8')],Response=stat1$Response[which(stat1$CellType=='CD8')],
  INVAX=stat1$Patient[which(stat1$CellType=='CD8')])
df_ratio=df_ratio[which(!is.na(df_ratio$Response)),]
df_ratio$Response=as.character(df_ratio$Response)
df_ratio=df_ratio[which(df_ratio$INVAX!='33'),]

ggpaired(df_ratio,x='Timepoint',y='CD8toTreg',fill = 'Timepoint', line.color = "gray", line.size = 0.4,id='INVAX')+
  scale_fill_manual(values=c( "#DDCC77","#AA4499"))+xlab(label='')+facet_wrap(.~Response)+
  stat_compare_means(paired = TRUE,label.y = 15,label = 'p.format')+ylab('CD8 to Treg ratio')+theme_classic(base_size = 16)+
  theme(legend.position = 'none')
compare_means(CD8toTreg~Response,df_ratio,group.by ='Timepoint')

