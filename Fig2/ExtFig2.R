#############################################Figure2c-e boxplots comparing cell type frequencies between response group#############

library(readxl)
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)


stat1=read_xlsx('~/Documents/ChenLab/Projects/HeadNeck/final.figures/Manuscript/1207/Supplementary Tables 6-13_v3.xlsx',
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
###########################Extended Figure 2#################################
#comparing between response at baseline
df_plot=stat1[which(!is.na(stat1$Response)&!is.na(stat1$Proportion)&stat1$Timepoint=='Baseline'),]
stat.test=compare_means(Proportion~Response,df_plot,p.adjust.method = 'fdr',group.by ='CellType')
stat.test=stat.test %>% add_xy_position(data = df_plot,formula = Proportion ~ Response)
ggboxplot(df_plot,x='Response',y='Proportion',fill='Response',add='point')+
  facet_grid(.~CellType,scales = 'free')+scale_fill_manual(values=c( "darkblue","skyblue"))+xlab(label='')+
  stat_pvalue_manual(stat.test, label = "p.adj")+
  theme_classic()+ ggtitle('Baseline')+
  theme(plot.title = element_text(hjust = 0.5),legend.position = 'none')
stat.test.baseline=stat.test

#comparing between response on induction
df_plot=stat1[which(!is.na(stat1$Response)&!is.na(stat1$Proportion)&stat1$Timepoint=='OnInduction'),]
stat.test=compare_means(Proportion~Response,df_plot,p.adjust.method = 'fdr',group.by ='CellType')
stat.test=stat.test %>% add_xy_position(data = df_plot,formula = Proportion ~ Response)
ggboxplot(df_plot,x='Response',y='Proportion',fill='Response',add='point')+
  facet_grid(.~CellType,scales = 'free')+scale_fill_manual(values=c( "darkblue","skyblue"))+xlab(label='')+stat_pvalue_manual(stat.test, label = "p.adj")+
  theme_classic()+ggtitle('OnInduction')+
  theme(plot.title = element_text(hjust = 0.5),legend.position = 'none')
stat.test.oi=stat.test

#compare between timepoint for MHR
df_paired=stat1[which(!is.na(stat1$Proportion)&stat1$Response=='MHR'),]
df_paired=df_paired[which(df_paired$Patient!='33'),]#patient 33 does not have baseline sample
stat.test=compare_means(Proportion ~ Timepoint,df_paired,p.adjust.method = 'fdr',group.by =c('CellType','Response'),paired = TRUE)
stat.test=stat.test %>% add_xy_position(data = df_paired,formula = Proportion ~ Timepoint)
ggpaired(df_paired,x='Timepoint',y='Proportion',fill = 'Timepoint', line.color = "gray", line.size = 0.4,id='Patient')+ggtitle('MHR')+
  facet_grid(Response~CellType,scales='free')+scale_fill_manual(values=c( "#DDCC77","#AA4499"))+xlab(label='')+stat_pvalue_manual(stat.test, label = "p.adj")+
  theme(plot.title = element_text(hjust = 0.5),legend.position = 'none')
stat.test.MHR=stat.test


#compare between timepoint for nonMHR
df_paired=stat1[which(!is.na(stat1$Proportion)&stat1$Response=='nMHR'),]
stat.test=compare_means(Proportion ~ Timepoint,df_paired,p.adjust.method = 'fdr',group.by =c('CellType','Response'),paired = TRUE)
stat.test=stat.test %>% add_xy_position(data = df_paired,formula = Proportion ~ Timepoint)
ggpaired(df_paired,x='Timepoint',y='Proportion',fill = 'Timepoint', line.color = "gray", line.size = 0.4,id='Patient')+ggtitle('MHR')+
  facet_grid(Response~CellType,scales='free')+scale_fill_manual(values=c( "#DDCC77","#AA4499"))+xlab(label='')+stat_pvalue_manual(stat.test, label = "p.adj")+
  theme(plot.title = element_text(hjust = 0.5),legend.position = 'none')
stat.test.nMHR=stat.test

