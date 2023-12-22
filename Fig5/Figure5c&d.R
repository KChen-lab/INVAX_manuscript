##############################Figure 5c&d##########################

library(ggplot2)
library(rstatix)
library(ggpubr)

#This code needs:
#scTCR-seq dataset includes cell barcode pseudoclonotype,clonotype (defined by CDR3A+CDR3B) and 
#whether clonotypes has expanded or not; the expansion test is done in Clonotyope_Expansion.R 
#cd8/cd4 meta data includes patient id, response,timepoint


meta=readRDS('./CD8_metadata.rds')
meta.cd4=readRDS('./CD4_metadata.rds')
meta_total=rbind(meta.cd4[,intersect(colnames(meta.cd4),colnames(meta))],meta[,intersect(colnames(meta.cd4),colnames(meta))])
wholevdj.dat=readRDS('./VDJ_data.rds')


meta=meta_cd8[,c('meta.INVAX.Study.ID','Timepoint','PATH.Response.new')]
meta_total=meta_t[,c('meta.INVAX.Study.ID','Timepoint','PATH.Response.new')]
colnames(wholevdj.dat)
ind1=match(intersect(wholevdj.dat$barcode,rownames(meta_total)),wholevdj.dat$barcode)
ind2=match(intersect(wholevdj.dat$barcode,rownames(meta_total)),rownames(meta_total))
identical(wholevdj.dat$barcode[ind1],rownames(meta_total)[ind2])
meta_total$Expansion_annot=NA
meta_total[ind2,'Expansion_annot']=wholevdj.dat[ind1,'Expansion_annot']
meta_total[ind2,'pseudoclonotype']=wholevdj.dat$pseudo_clonotype[ind1]
meta_total[ind2,'vdj_barcode']=wholevdj.dat$barcode[ind1]

#map scTCR data to cd8 data
ind1=match(intersect(wholevdj.dat$barcode,rownames(meta)),wholevdj.dat$barcode)
ind2=match(intersect(wholevdj.dat$barcode,rownames(meta)),rownames(meta))
identical(wholevdj.dat$barcode[ind1],rownames(meta)[ind2])
meta$Expansion_annot=NA
meta[ind2,'Expansion_annot']=wholevdj.dat[ind1,'Expansion_annot']
meta[ind2,'pseudoclonotype']=wholevdj.dat$pseudo_clonotype[ind1]
meta[ind2,'vdj_barcode']=wholevdj.dat$barcode[ind1]

pts=unique(meta$patient)
df=data.frame(invax=pts)
for(i in 1:nrow(df)){
  patient_data=meta[which(meta$patient==df$invax[i]),]
  patient_data_all=meta_total[which(meta_total$patient==df$invax[i]),]
  df[i,'Response']=unique(meta$Response[which(meta$patient==df$invax[i])])
  df[i,'expanded_count']=length(unique(patient_data$pseudoclonotype[which(patient_data$Expansion_annot=='Expanded_Sig')]))
  df[i,'contracted_count']=length(unique(patient_data$pseudoclonotype[which(patient_data$Expansion_annot=='Contracted_Sig')]))
  df[i,'novel_count']=length(unique(patient_data$pseudoclonotype[which(patient_data$Expansion_annot=='Novel_Sig')]))
  df[i,'expanded_fraction']=length(unique(patient_data$pseudoclonotype[which(patient_data$Expansion_annot=='Expanded_Sig')]))/length(unique(patient_data_all$pseudoclonotype[which(patient_data_all$Timepoint==1)]))*100
  df[i,'contracted_fraction']=length(unique(patient_data$pseudoclonotype[which(patient_data$Expansion_annot=='Contracted_Sig')]))/length(unique(patient_data_all$pseudoclonotype[which(patient_data_all$Timepoint==1)]))*100
  df[i,'novel_fraction']=length(unique(patient_data$pseudoclonotype[which(patient_data$Expansion_annot=='Novel_Sig')]))/length(unique(patient_data_all$pseudoclonotype[which(patient_data_all$Timepoint==2)]))*100
}
df=df[which(df$invax!='INVAX026'&df$invax!='INVAX033'&!is.na(df$Response)),]

p1=ggplot(df,aes(x=Response,y=expanded_fraction,fill=Response))+geom_boxplot()+geom_point()+
  stat_compare_means()+scale_fill_manual(values=c('royalblue','lightblue'))+theme_classic(base_size = 15)+
  theme(legend.position = 'none',axis.text=element_text(color='black'))+ylab('% Fraction of expanded clones')+xlab('')
p2=ggplot(df,aes(x=Response,y=contracted_fraction,fill=Response))+geom_boxplot()+geom_point()+
  stat_compare_means()+scale_fill_manual(values=c('royalblue','lightblue'))+theme_classic(base_size = 15)+
  theme(legend.position = 'none',axis.text=element_text(color='black'))+ylab('% Fraction of contracted clones')+xlab('')
# p3=ggplot(df,aes(x=Response,y=novel_fraction,fill=Response))+geom_boxplot()+geom_point()+
#   stat_compare_means()+scale_fill_manual(values=c('royalblue','lightblue'))+theme_classic(base_size = 15)+
#   theme(legend.position = 'none',axis.text=element_text(color='black'))+ylab('% Fraction of novel clones')+xlab('')
# p1+p2+p3
p1+p2

stat.test=rbind(wilcox_test(expanded_fraction~Response,data=df),
                wilcox_test(contracted_fraction~Response,data=df))

