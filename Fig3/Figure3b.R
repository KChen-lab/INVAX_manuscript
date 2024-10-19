############################################# Figure3a boxplots comparing effector Treg cell subtype frequencies between response group#############

library(readxl)
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)


stat1=read_xlsx('./Supplementary Tables 6-13_v3.xlsx',
                sheet='Supplementary Table 7',range = 'A3:AE73')
response=read_xlsx('../DataShare/response.xlsx')
response=as.data.frame(response)

stat1[,3:31]=100*stat1[,3:31]/rowSums(stat1[,3:31])

stat1=melt(stat1,id.vars = c('Patient ID','Timepoint'),variable.name ='CellSubset')
colnames(stat1)=c('Patient','Timepoint','CellSubset','Proportion')
stat1$Response=as.character(response[match(stat1$Patient,response$Patient),'Response'])
stat1$Response=factor(stat1$Response,levels=c('MHR','nMHR'))
stat1=stat1[which(!is.na(stat1$Response)),]
stat1=stat1[which(stat1$Patient!='26'),]
stat1=stat1[which(stat1$Patient!='33'),]

###########################Extended Figure 3c-d#################################
#comparing between response at baseline
df_plot=stat1[which(!is.na(stat1$Response)&!is.na(stat1$Proportion)&stat1$Timepoint=='Baseline'),]
cells=as.character(unique(df_plot$CellSubset))[1:14]
df_plot=df_plot[which(!is.na(match(df_plot$CellSubset,cells))),]
df_plot$CellSubset=factor(df_plot$CellSubset,levels = cells)
stat.test=compare_means(Proportion~Response,df_plot,p.adjust.method = 'fdr',group.by ='CellSubset')
stat.test=stat.test %>% add_xy_position(data = df_plot,formula = Proportion ~ Response)
stat.test.baseline=stat.test

#comparing between response on induction
df_plot=stat1[which(!is.na(stat1$Response)&!is.na(stat1$Proportion)&stat1$Timepoint=='OnInduction'),]
cells=as.character(unique(df_plot$CellSubset))[1:14]
df_plot=df_plot[which(!is.na(match(df_plot$CellSubset,cells))),]
df_plot$CellSubset=factor(df_plot$CellSubset,levels = cells)
stat.test=compare_means(Proportion~Response,df_plot,p.adjust.method = 'fdr',group.by ='CellSubset')
stat.test=stat.test %>% add_xy_position(data = df_plot,formula = Proportion ~ Response)
stat.test.oi=stat.test

#compare between timepoints for MHR-paired
df_paired=stat1[which(!is.na(stat1$Proportion)&stat1$Response=='MHR'),]
cells=as.character(unique(df_paired$CellSubset))[1:14]
df_paired=df_paired[which(!is.na(match(df_paired$CellSubset,cells))),]
stat.test=compare_means(Proportion ~ Timepoint,df_paired,p.adjust.method = 'fdr',group.by =c('CellSubset','Response'),paired = TRUE)
stat.test=stat.test %>% add_xy_position(data = df_paired,formula = Proportion ~ Timepoint)
stat.test.MHR=stat.test


#compare between timepoints for nonMHR-paired
df_paired=stat1[which(!is.na(stat1$Proportion)&stat1$Response=='nMHR'),]
cells=as.character(unique(df_paired$CellSubset))[1:14]
df_paired=df_paired[which(!is.na(match(df_paired$CellSubset,cells))),]
stat.test=compare_means(Proportion ~ Timepoint,df_paired,p.adjust.method = 'fdr',group.by =c('CellSubset','Response'),paired = TRUE)
stat.test=stat.test %>% add_xy_position(data = df_paired,formula = Proportion ~ Timepoint)
stat.test.nMHR=stat.test



plot_celltype_proportion<-function(data,celltype){
  df_malig=data[which(!is.na(data$Proportion)&data$CellSubset==celltype&!is.na(data$Response)),]
  pvalues <- data.frame(
    Response = c("MHR", "nMHR",'MHR','nMHR'),  
    x = c(1.5,1.5,2,1.5),  
    y = c(round(max(df_malig$Proportion))-5,round(max(df_malig$Proportion))-5,round(max(df_malig$Proportion))+5,round(max(df_malig$Proportion))+5), 
    label = c(formatC(stat.test.MHR$p.adj[which(stat.test.MHR$CellSubset==celltype)], digits=2),
              formatC(stat.test.nMHR$p.adj[which(stat.test.nMHR$CellSubset==celltype)], digits=2),
              paste0('B:',formatC(stat.test.baseline$p.adj[which(stat.test.baseline$CellSubset==celltype)], digits=2)),
              paste0('OI:',formatC(stat.test.oi$p.adj[which(stat.test.oi$CellSubset==celltype)], digits=2))
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

p=plot_celltype_proportion(stat1,'Treg.Effector')
p
