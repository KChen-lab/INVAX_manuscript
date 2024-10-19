
###########################################Extended Figure 7b

library(stringr)
library(ggplot2)
library(rstatix)
library(ggpubr)


#using proportion_vdj
get_percent_occupancy <- function(patient_data,timepoint,percent){
  patient_data_pre=patient_data[which(patient_data$Timepoint==timepoint),]
  tt=table(patient_data_pre$pseudoclonotype,patient_data_pre$proportion_cdr3.aa_local)
  tt=as.data.frame(tt)
  tt=tt[which(tt$Freq!=0),]
  if(nrow(tt)>0){
    tt=tt[order(tt$Var2,decreasing = T),]
    count=0
    for( i in 1:nrow(tt)){
      count=count+as.numeric(as.character(tt$Var2[i]))
      if (count>0.1){
        no_clones=i
        break
      }
    }
  }
  if(nrow(tt)==0){
    no_clones=NA
  }
  return(no_clones)
}


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



percent=0.1
df=data.frame(invax=unique(meta$patient))
tt=as.data.frame(table(meta$patient,meta$Response))
tt=tt[which(tt$Freq!=0),]
rownames(tt)=tt$Var1
df$response=tt[df$invax,'Var2']

for(i in 1:nrow(df)){
  df[i,'pre_count']=get_percent_occupancy(meta[which(meta$patient==df$invax[i]&!is.na(meta$vdj_barcode)),],1,percent)

}
df1=df[which(!is.na(df$response)&df$invax!='INVAX026'&df$invax!='INVAX033'),]
ggplot(df1,aes(x=response,y=pre_count,fill=response))+geom_boxplot() +
  geom_point()+stat_compare_means()+theme_classic(base_size = 15)+
  scale_fill_manual(values=c("darkblue","skyblue"))+ylab('No of clones occupaying \n 10% of whole basal TCR repertoire')+
  xlab('')+theme(legend.position = 'none')

stat.test=wilcox_test(pre_count~response,data=df1)
