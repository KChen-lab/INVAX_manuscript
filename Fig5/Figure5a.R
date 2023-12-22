
###########################################Figure 5a
library(ggplot2)
library(ggpubr)
library(rstatix)

get_percent<- function(patient_data,timepoint,percent){
  patient_data_pre=patient_data[which(patient_data$Timepoint==timepoint),]
  if(nrow(patient_data_pre)>0){ 
    tt=as.data.frame(table(patient_data_pre$pseudoclonotype,patient_data_pre$proportion_cdr3.aa_local))
    tt=tt[which(tt$Freq!=0),]
    no_clones=length(which(as.numeric(as.character(tt$Var2)) > percent|as.numeric(as.character(tt$Var2)) == percent))
  }
  if(nrow(patient_data_pre)==0){
    no_clones=NA
  }
  return(no_clones)
}



#This code needs:
#scTCR-seq dataset includes cell barcode pseudoclonotype,clonotype (defined by CDR3A+CDR3B) and frequency 
#cd8/cd4 meta data includes patient id, response,timepoint

wholevdj.dat=readRDS('./VDJ_data.rds')
meta=readRDS('./CD8_metadata.rds')
meta.cd4=readRDS('./CD4_metadata.rds')
meta_total=rbind(meta.cd4[,intersect(colnames(meta.cd4),colnames(meta))],meta[,intersect(colnames(meta.cd4),colnames(meta))])

temp=intersect(wholevdj.dat$barcode,rownames(meta_total))
ind1=match(temp,wholevdj.dat$barcode)
ind2=match(temp,rownames(meta_total))
identical(wholevdj.dat$barcode[ind1],rownames(meta_total)[ind2])
meta_total[ind2,'pseudoclonotype']=wholevdj.dat$pseudo_clonotype[ind1]
meta_total[ind2,'proportion_cdr3.aa_local']=wholevdj.dat$proportion_cdr3.aa_local[ind1]
meta_total[ind2,'vdj_barcode']=wholevdj.dat$barcode[ind1]

temp=intersect(wholevdj.dat$barcode,rownames(meta))
ind1=match(temp,wholevdj.dat$barcode)
ind2=match(temp,rownames(meta))
identical(wholevdj.dat$barcode[ind1],rownames(meta)[ind2])
meta[ind2,'pseudoclonotype']=wholevdj.dat$pseudo_clonotype[ind1]
meta[ind2,'proportion_cdr3.aa_local']=wholevdj.dat$proportion_cdr3.aa_local[ind1]
meta[ind2,'vdj_barcode']=wholevdj.dat$barcode[ind1]

percent=0.01
df=data.frame(invax=unique(meta$patient))
tt=as.data.frame(table(meta$patient,meta$Response))
tt=tt[which(tt$Freq!=0),]
rownames(tt)=tt$Var1
df$response=tt[df$invax,'Var2']
for(i in 1:nrow(df)){
  df[i,'pre_count']=get_percent(meta[which(meta$patient==df$invax[i]&!is.na(meta$vdj_barcode)),],1,percent)
  df[i,'post_count']=get_percent(meta[which(meta$patient==df$invax[i]&!is.na(meta$vdj_barcode)),],2,percent)
  df[i,'pre_count_total']=length(unique(meta_total$pseudoclonotype[which(meta_total$patient==df$invax[i]&meta_total$Timepoint==1&!is.na(meta_total$vdj_barcode))]))
  df[i,'post_count_total']=length(unique(meta_total$pseudoclonotype[which(meta_total$patient==df$invax[i]&meta_total$Timepoint==2&!is.na(meta_total$vdj_barcode))]))
}
df1=df[which(!is.na(df$response)),]
df1$pre_percent=df1$pre_count/df1$pre_count_total*100
df1$post_percent=df1$post_count/df1$post_count_total*100

df1=df1[which(df1$invax!='INVAX026'&df1$invax!='INVAX033'),] #BECASUE INVAX026 only has 13 t cells at baseline

ggplot(df1,aes(x=response,y=pre_percent,fill=response))+geom_boxplot() +
  geom_point()+stat_compare_means()+theme_classic(base_size = 15)+
  scale_fill_manual(values=c("royalblue","lightblue"))+ylab('Proportion of unique clones > 1%\n at baseline')+
  xlab('')+theme(legend.position = 'none')

stat.test=wilcox_test(pre_percent~response,data=df1)
