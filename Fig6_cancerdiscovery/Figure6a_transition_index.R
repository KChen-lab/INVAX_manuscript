library(philentropy)
library(ggplot2)
library(ggpubr)
library(rstatix)


######################################################Figure 6a  Transition index between MHR and nMHR###########################
#This code needs dataset file including patient id ('meta.INVAX.Study.ID'), Timepoint, cell subset annotation, 
#clonotype id or cdr3a+cdr3b sequence, whether patient have paired timepoints ('Match') for each CD8+ T cell.
#All of the information and data can be download through XXXXX.


meta=readRDS('./CD8_metadata.rds')
wholevdj.dat=readRDS('./VDJ_data.rds')

temp=intersect(wholevdj.dat$barcode,rownames(meta))
ind1=match(temp,wholevdj.dat$barcode)
ind2=match(temp,rownames(meta))
identical(wholevdj.dat$barcode[ind1],rownames(meta)[ind2])
meta[ind2,'clonotype']=wholevdj.dat$pseudo_clonotype[ind1]
meta[ind2,'vdj_barcode']=wholevdj.dat$barcode[ind1]

colnames(meta)
meta=meta[which(!is.na(meta$vdj_barcode)),]#remove T cells without VDJ information

tt=table(meta$patient,meta$Timepoint)
pts=names(which(tt[,1]!=0&tt[,2]))

df=data.frame(patient=pts)
for(i in 1:nrow(df)){
  kl=0
  temp=meta[which(meta$patient==pts[i]),]
  clones=unique(temp$clonotype)
  for(j in 1:length(clones)){
    temp.1=temp[which(temp$clonotype==clones[j]),]
    temp.2=table(temp.1$Timepoint)
   length(temp.2)
    if(length(temp.2)>1){
      temp.3=table(temp.1$Timepoint,temp.1$CellSubset)
      temp_pre=temp.3[1,]
      temp_post=temp.3[2,]
      kl.temp=KL(rbind(temp_pre,temp_post),est.prob = 'empirical')
      p.temp=temp.2['1']/length(which(temp$Timepoint==1))
      kl=kl+kl.temp*p.temp
    }
  }
    df[i,'Transition']=kl
    df[i,'Response']=unique(temp$Response)
}

df_plot=df[which(!is.na(df$Response)),]
ggplot(df_plot,aes(x=Response,y=Transition,fill=Response))+geom_boxplot()+
  stat_compare_means()+geom_point()+theme_bw()+scale_fill_manual(values=c('royalblue','lightblue'))

stats=wilcox_test(Transition~Response, data=df_plot)
