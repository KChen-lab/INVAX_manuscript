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

########################### Figure 2c-e#################################
#comparing between response at baseline
df_plot=stat1[which(!is.na(stat1$Response)&!is.na(stat1$Proportion)&stat1$Timepoint=='Baseline'),]
stat.test=compare_means(Proportion~Response,df_plot,p.adjust.method = 'fdr',group.by ='CellType')
stat.test.baseline=stat.test

#comparing between response on induction
df_plot=stat1[which(!is.na(stat1$Response)&!is.na(stat1$Proportion)&stat1$Timepoint=='OnInduction'),]
stat.test=compare_means(Proportion~Response,df_plot,p.adjust.method = 'fdr',group.by ='CellType')
stat.test=stat.test %>% add_xy_position(data = df_plot,formula = Proportion ~ Response)
stat.test.oi=stat.test

#compare between timepoint for MHR
df_paired=stat1[which(!is.na(stat1$Proportion)&stat1$Response=='MHR'),]
df_paired=df_paired[which(df_paired$Patient!='33'),]#patient 33 does not have baseline sample
stat.test=compare_means(Proportion ~ Timepoint,df_paired,p.adjust.method = 'fdr',group.by =c('CellType','Response'),paired = TRUE)
stat.test=stat.test %>% add_xy_position(data = df_paired,formula = Proportion ~ Timepoint)
stat.test.MHR=stat.test


#compare between timepoint for nonMHR
df_paired=stat1[which(!is.na(stat1$Proportion)&stat1$Response=='nMHR'),]
stat.test=compare_means(Proportion ~ Timepoint,df_paired,p.adjust.method = 'fdr',group.by =c('CellType','Response'),paired = TRUE)
stat.test=stat.test %>% add_xy_position(data = df_paired,formula = Proportion ~ Timepoint)
stat.test.nMHR=stat.test



plot_celltype_proportion<-function(data,celltype){
  df_malig=data[which(!is.na(data$Proportion)&data$CellType==celltype&!is.na(data$Response)&data$Patient!=33),]
  pvalues <- data.frame(
    Response = c("MHR", "nMHR",'MHR','nMHR'),  
    x = c(1.5,1.5,2,1.5),  
    y = c(round(max(df_malig$Proportion))-5,round(max(df_malig$Proportion))-5,round(max(df_malig$Proportion))+5,round(max(df_malig$Proportion))+5), 
    label = c(formatC(stat.test.MHR$p.adj[which(stat.test.MHR==celltype)], digits=2),
              formatC(stat.test.nMHR$p.adj[which(stat.test.nMHR==celltype)], digits=2),
              paste0('B:',formatC(stat.test.baseline$p.adj[which(stat.test.baseline==celltype)], digits=2)),
              paste0('OI:',formatC(stat.test.oi$p.adj[which(stat.test.oi==celltype)], digits=2))
    )  
  )
  p=ggpaired(df_malig,x='Timepoint',y='Proportion',fill = 'Timepoint', line.color = "gray", line.size = 0.4,id='Patient')+
    scale_fill_manual(values=c( "#DDCC77","#AA4499"))+xlab(label='')+
    theme(plot.title = element_text(hjust = 0.5),legend.position = 'none')+
    facet_wrap(Response~.)+
    geom_text(data = pvalues, aes(x = x, y = y, label = label), size = 4, inherit.aes = FALSE)+
    ggtitle(celltype)+ylab('%Total TIME cells')
  return(p)
}

p.cd8=plot_celltype_proportion(stat1,'CD8')
p.cd8
ggsave('~/Documents/ChenLab/Projects/HeadNeck/final.figures/Figure2/CD8.pdf',height = 4,width = 5)
p.cd4=plot_celltype_proportion(stat1,'CD4')
p.cd4
ggsave('~/Documents/ChenLab/Projects/HeadNeck/final.figures/Figure2/CD4.pdf',height = 4,width = 5)

p.malig=plot_celltype_proportion(stat1,'Malignant')
p.malig  
ggsave('~/Documents/ChenLab/Projects/HeadNeck/final.figures/Figure2/Malignant.pdf',height = 4,width = 5)
