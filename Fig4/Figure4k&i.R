
library(readxl)
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(reshape2)

stat1=read_xlsx('~/Documents/ChenLab/Projects/HeadNeck/final.figures/Manuscript/1207/Supplementary Tables 6-13_v3.xlsx',
                sheet='Supplementary Table 7',range = 'A3:AE73')
stat1=as.data.frame(stat1)
response=read_xlsx('~/Documents/ChenLab/Projects/HeadNeck/final.figures/Manuscript/Data and code/DataShare/response.xlsx')
response=as.data.frame(response)

stat1[,3:31]=stat1[,3:31]/rowSums(stat1[,3:31])
stat1$Response=as.character(response[match(stat1$`Patient ID`,response$Patient),'Response'])
stat1$Response=factor(stat1$Response,levels=c('MHR','nMHR'))
colnames(stat1)

nonpaired=stat1$`Patient ID`[which(is.na(stat1$CD4.Naive))]
stat1=stat1[which(is.na(match(stat1$`Patient ID`,nonpaired))),] #remove non paired sample
stat1=stat1[which(stat1$`Patient ID`!=26),] #remove patient 26

pre=stat1[which(stat1$Timepoint=='Baseline'),]
post=stat1[which(stat1$Timepoint=='OnInduction'),]

#all CD8 clusters
cells=colnames(stat1)[17:31]
method='spearman'
df=data.frame(celltype=cells)
for(i in 1:15){
  temp=cor.test((post[,cells[i]]-pre[,cells[i]])/pre[,cells[i]],
                (post[,'Treg.Effector']-pre[,'Treg.Effector'])/pre[,'Treg.Effector'],method=method)
  df[i,'delta_spearman_cor']=temp$estimate
  df[i,'delta_spearman_pval']=temp$p.value
  df[i,'delta_spearman_statistic']=temp$statistic
}
colnames(df)

df_plot=data.frame(celltype=df$celltype,
                   cor=df[,'delta_spearman_cor'],
                   raw_pvalue=as.numeric(format(df[,'delta_spearman_pval'],digits = 4,scientific = F)),
                   pvals=as.numeric(format(p.adjust(df[,'delta_spearman_pval'],method='fdr'),digits = 4,scientific = T)),
                   statistic=df$delta_spearman_statistic)

df_plot=df_plot[order(df_plot$cor),]
df_plot$celltype=factor(df_plot$celltype,levels=rev(df_plot$celltype))
df_plot$pval.postition=1
df_plot$pval.postition[which(df_plot$cor<0)]=0.1
df_plot$pval.postition[which(df_plot$cor>0)]=-0.1
df_plot$pval.symbol=''
df_plot$pval.symbol[which(df_plot$pvals<0.05)]='*'
df_plot$pval.symbol[which(df_plot$pvals<0.01)]='**'
df_plot$pval.symbol[which(df_plot$pvals<0.001)]='***'
ggplot(df_plot[which(df_plot$celltype!='CD4CD8'&df_plot$celltype!='CD8.Undefined'),],aes(x=cor,y=celltype,label=pval.symbol))+geom_bar(stat = 'identity',fill='grey70',color='black')+theme_bw()+xlab('Spearman Rho')+
  geom_text(aes(x=pval.postition,y=celltype),size=10)+ylab('')+
  theme(axis.text=element_text(color='black'))

df_plot2=data.frame(Treg=(post$Treg.Effector-pre$Treg.Effector)*100,IFNG=(post$`CD8.Effector.1(IFNG+)`-pre$`CD8.Effector.1(IFNG+)`)*100)
rownames(df_plot2)=rownames(post)
df_plot2
ggplot(df_plot2,aes(y=IFNG,x=Treg))+geom_point(size=3)+theme_classic(base_size = 15)+
  geom_smooth(method='lm')+ggtitle(paste0('Spearman Rho: ',format(df$delta_spearman_cor[4],digits = 4)))+xlab('Delta proportion of \n Treg.Effector')+
  ylab('Delta proportion of \n CD8.Effector.1(IFNG+)')+theme(axis.text=element_text(color='black'))
ggsave('/Volumes/HTS-Intel-4/Research/HeadandNeck/hn.paper/final.figures/Figure4/Scatter_plot_treg_ifng_deta.pdf',width = 5.5,height = 5)



