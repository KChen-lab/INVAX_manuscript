#############################################Figure2b scatter plot for tumor viable change and malignant cell frequency change#############

library(ggbreak)
library(readxl)
library(ggplot2)


stat1=read_xlsx('~/Documents/ChenLab/Projects/HeadNeck/final.figures/Manuscript/1207/Supplementary Tables 6-13_v3.xlsx',
               sheet='Supplementary Table 6',range = 'A3:N73')
response=read_xlsx('../DataShare/response.xlsx')
response=as.data.frame(response)

stat1[,4:14]=stat1[,4:14]/stat1$`Total No. Cells after QC`
stat1=stat1[,-3]
colnames(stat1)
stat1=melt(stat1,id.vars = c('Patient ID','Timepoint'),variable.name ='CellType')
colnames(stat1)=c('Patient','Timepoint','CellType','Proportion')

df_malig=stat1[which(stat1$CellType=='Malignant'&stat1$Timepoint=='Baseline'),]
colnames(df_malig)[4]='Pre'
df_malig=df_malig[,-c(2,3)]
df_malig$Post=stat1$Proportion[which(stat1$CellType=='Malignant'&stat1$Timepoint=='OnTreatment')]

df_malig$ViableChange=as.numeric(as.character(response[match(df_malig$Patient,response$Patient),'ViableChange']))
df_malig$Response=as.character(response[match(df_malig$Patient,response$Patient),'Response'])

df_malig$Reduction=df_malig$Post-df_malig$Pre
df_malig$Reduction_normalized=(df_malig$Post-df_malig$Pre)/df_malig$Pre*100

df_malig=df_malig[which(df_malig$Patient!=33),]
#write_xlsx(df_malig,'/Volumes/HTS-Intel-4/Research/HeadandNeck/hn.paper/final.figures/Figure2/SourceData/SourceData_Figure2C-scatter_viablechange_malignant_reduction_spearman.xlsx')
test=cor.test(df_malig$ViableChange[which(!is.na(df_malig$ViableChange))],df_malig$Reduction_normalized[which(!is.na(df_malig$ViableChange))],
         method = 'spearman')
test$p.value
test$estimate

ggplot(df_malig[which(!is.na(df_malig$ViableChange)),],aes(x=Reduction_normalized,y=ViableChange,color=Response))+geom_point(size=5,alpha=0.5)+
  theme_classic(base_size = 16)+xlab('Single cell malignant cell reduction (Normalized) %')+
  ylab('Tumor viability change (Normalized) %')+scale_x_break(c(350,1500),ticklabels = c(1500,1600))+
  scale_color_manual(values=c('darkblue','skyblue'))+theme(plot.title = element_text(hjust = 0.5))



