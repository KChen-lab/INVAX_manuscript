#############################################Extended Figure3c-d boxplots comparing CD4+ T cell subtype frequencies between response group#############

library(readxl)
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(reshape2)

stat1=read_xlsx('~/Documents/ChenLab/Projects/HeadNeck/final.figures/Manuscript/1207/Supplementary Tables 6-13_v3.xlsx',
                sheet='Supplementary Table 7',range = 'A3:AE73')
response=read_xlsx('../DataShare/response.xlsx')
response=as.data.frame(response)

stat1[,3:31]=100*stat1[,3:31]/rowSums(stat1[,3:31])

stat1=melt(stat1,id.vars = c('Patient ID','Timepoint'),variable.name ='CellSubset')
colnames(stat1)=c('Patient','Timepoint','CellSubset','Proportion')
stat1$Response=as.character(response[match(stat1$Patient,response$Patient),'Response'])
stat1$Response=factor(stat1$Response,levels=c('MHR','nMHR'))

###########################Extended Figure 3c-d#################################
#comparing between response at baseline
df_plot=stat1[which(!is.na(stat1$Response)&!is.na(stat1$Proportion)&stat1$Timepoint=='Baseline'),]
df_plot=df_plot[which(df_plot$Patient!=26),]#patient 26 is removed due to small number of total T cells at baseline
cells=as.character(unique(df_plot$CellSubset))[1:14]
df_plot=df_plot[which(!is.na(match(df_plot$CellSubset,cells))),]
df_plot$CellSubset=factor(df_plot$CellSubset,levels = cells)
stat.test=compare_means(Proportion~Response,df_plot,p.adjust.method = 'fdr',group.by ='CellSubset')
stat.test=stat.test %>% add_xy_position(data = df_plot,formula = Proportion ~ Response)

ggboxplot(df_plot,x='Response',y='Proportion',fill='Response',add='point')+
  facet_grid(.~CellSubset,scales = 'free')+scale_fill_manual(values=c( "royalblue","skyblue"))+xlab(label='')+
  stat_pvalue_manual(stat.test, label = "p.adj")+
  theme_classic()+ ggtitle('Baseline')+
  theme(plot.title = element_text(hjust = 0.5),legend.position = 'none')
stat.test.baseline=stat.test

#comparing between response on induction
df_plot=stat1[which(!is.na(stat1$Response)&!is.na(stat1$Proportion)&stat1$Timepoint=='OnInduction'),]
cells=as.character(unique(df_plot$CellSubset))[1:14]
df_plot=df_plot[which(!is.na(match(df_plot$CellSubset,cells))),]
df_plot$CellSubset=factor(df_plot$CellSubset,levels = cells)
stat.test=compare_means(Proportion~Response,df_plot,p.adjust.method = 'fdr',group.by ='CellSubset')
stat.test=stat.test %>% add_xy_position(data = df_plot,formula = Proportion ~ Response)
ggboxplot(df_plot,x='Response',y='Proportion',fill='Response',add='point')+
  facet_grid(.~CellSubset,scales = 'free')+scale_fill_manual(values=c( "royalblue","skyblue"))+xlab(label='')+
  stat_pvalue_manual(stat.test[1:14,], label = "p.adj")+
  theme_classic()+ ggtitle('Baseline')+
  theme(plot.title = element_text(hjust = 0.5),legend.position = 'none')
stat.test.oi=stat.test

#compare between timepoints for MHR-paired
df_paired=stat1[which(!is.na(stat1$Proportion)&stat1$Response=='MHR'),]
df_paired=df_paired[which(df_paired$Patient!='33'),]#patient 33 does not have baseline sample
cells=as.character(unique(df_paired$CellSubset))[1:14]
df_paired=df_paired[which(!is.na(match(df_paired$CellSubset,cells))),]
stat.test=compare_means(Proportion ~ Timepoint,df_paired,p.adjust.method = 'fdr',group.by =c('CellSubset','Response'),paired = TRUE)
stat.test=stat.test %>% add_xy_position(data = df_paired,formula = Proportion ~ Timepoint)
ggpaired(df_paired,x='Timepoint',y='Proportion',fill = 'Timepoint', line.color = "gray", line.size = 0.4,id='Patient')+ggtitle('MHR')+
  facet_grid(Response~CellSubset,scales='free')+scale_fill_manual(values=c( "#DDCC77","#AA4499"))+xlab(label='')+stat_pvalue_manual(stat.test, label = "p.adj")+
  theme(plot.title = element_text(hjust = 0.5),legend.position = 'none')
stat.test.MHR=stat.test


#compare between timepoints for nonMHR-paired
df_paired=stat1[which(!is.na(stat1$Proportion)&stat1$Response=='nMHR'),]
df_paired=df_paired[which(df_paired$Patient!='26'),]#patient 33 does not have baseline sample at baseline
cells=as.character(unique(df_paired$CellSubset))[1:14]
df_paired=df_paired[which(!is.na(match(df_paired$CellSubset,cells))),]
stat.test=compare_means(Proportion ~ Timepoint,df_paired,p.adjust.method = 'fdr',group.by =c('CellSubset','Response'),paired = TRUE)
stat.test=stat.test %>% add_xy_position(data = df_paired,formula = Proportion ~ Timepoint)
ggpaired(df_paired,x='Timepoint',y='Proportion',fill = 'Timepoint', line.color = "gray", line.size = 0.4,id='Patient')+ggtitle('nMHR')+
  facet_grid(Response~CellSubset,scales='free')+scale_fill_manual(values=c( "#DDCC77","#AA4499"))+xlab(label='')+stat_pvalue_manual(stat.test, label = "p.adj")+
  theme(plot.title = element_text(hjust = 0.5),legend.position = 'none')
stat.test.nMHR=stat.test


