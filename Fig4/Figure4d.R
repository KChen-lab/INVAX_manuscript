
################################Figure4d - TRM metacluster###########################
library(readxl)
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(reshape2)

stat1=read_xlsx('./Supplementary Tables 6-13_v3.xlsx',
                sheet='Supplementary Table 7',range = 'A3:AE73')
response=read_xlsx('../DataShare/response.xlsx')
response=as.data.frame(response)

stat1[,3:31]=100*stat1[,3:31]/rowSums(stat1[,3:31])
stat1$Response=as.character(response[match(stat1$`Patient ID`,response$Patient),'Response'])
stat1$Response=factor(stat1$Response,levels=c('MHR','nMHR'))
colnames(stat1)

TRM=grep('TRM',colnames(stat1))
TEM=grep('TEM',colnames(stat1))
Effector=grep('Effector.',colnames(stat1))

rowSums(stat1[,3:31])
stat2=stat1[,c(1,2,32)]
stat2$Other=rowSums(stat1[,setdiff(3:31,c(TRM,TEM,Effector))])
stat2$TRM=rowSums(stat1[,TRM])
stat2$TEM=rowSums(stat1[TEM])
stat2$Effector=rowSums(stat1[,Effector])

colnames(stat2)
stat2=melt(stat2,id.vars = c('Patient ID','Timepoint','Response'),variable.name ='metacluster')
colnames(stat2)=c('Patient','Timepoint','Response','metacluster','Proportion')


#comparing between response at baseline

df_plot=stat2[which(!is.na(stat2$Response)&!is.na(stat2$Proportion)&stat2$Timepoint=='Baseline'),]
df_plot=df_plot[which(df_plot$Patient!=26),]
df_plot=df_plot[which(df_plot$metacluster!='Other'),]

metacluster_type='TRM'

ggplot(df_plot[which(df_plot$metacluster==metacluster_type),],aes(x=Response,y=Proportion,fill=Response))+geom_boxplot()+geom_point()+
  stat_compare_means()+scale_

ggplot(df_plot[which(df_plot$metacluster==metacluster_type),],aes(x=Response,y=Proportion,fill=Response))+geom_boxplot()+geom_point(size=2)+
  stat_compare_means()+scale_fill_manual(values=c('royalblue','lightblue'))+theme_classic(base_size = 15)+ylab('% Total T cells')+xlab('')+
  theme(legend.position = 'none',axis.text = element_text(color='black'))
wilcox_test(Proportion~Response,data=df_plot[which(df_plot$metacluster==metacluster_type),],conf.level  = 0.95)

#this code can be also used to generate supplementary figure h by changing metacluster_type to TEM or Effector
