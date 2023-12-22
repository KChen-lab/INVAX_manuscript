

library(ggplot2)
library(rstatix)
library(ggpubr)
library(Seurat)

##################################################################Figure 5B##########################

#This code needs:
#scTCR-seq dataset includes cell barcode pseudoclonotype,clonotype (defined by CDR3A+CDR3B) and frequency 
#cd8 meta data includes patient id, response,timepoint

data=readRDS('./CD8.rds')
wholevdj.dat=readRDS('./VDJ_data.rds')
#{r clonotype size umap }
meta=data@meta.data

#map scTCRseq data to CD8
temp=intersect(wholevdj.dat$barcode,rownames(meta))
ind1=match(temp,wholevdj.dat$barcode)
ind2=match(temp,rownames(meta))
identical(wholevdj.dat$barcode[ind1],rownames(meta)[ind2])
meta[ind2,'pseudoclonotype']=wholevdj.dat$pseudo_clonotype[ind1]
meta[ind2,'proportion_cdr3.aa_local']=wholevdj.dat$proportion_cdr3.aa_local[ind1]
meta[ind2,'vdj_barcode']=wholevdj.dat$barcode[ind1]

meta$clone_size_cdr3_group_per_sample=NA
meta$clone_size_cdr3_group_per_sample[which(meta$proportion_cdr3.aa_local<0.0001)]='<0.01'
meta$clone_size_cdr3_group_per_sample[which(meta$proportion_cdr3.aa_local<0.001&(meta$proportion_cdr3.aa_local>0.0001|meta$proportion_cdr3.aa_local==0.0001))]='[0.01,0.1)'
meta$clone_size_cdr3_group_per_sample[which(meta$proportion_cdr3.aa_local<0.01&(meta$proportion_cdr3.aa_local>0.001|meta$proportion_cdr3.aa_local==0.001))]='[0.1,1)'
meta$clone_size_cdr3_group_per_sample[which(meta$proportion_cdr3.aa_local<0.1&(meta$proportion_cdr3.aa_local>0.01|meta$proportion_cdr3.aa_local==0.01))]='[1,10)'
meta$clone_size_cdr3_group_per_sample[which(meta$proportion_cdr3.aa_local>0.1|meta$proportion_cdr3.aa_local==0.1)]='>10'
meta$clone_size_cdr3_group_per_sample=factor(meta$clone_size_cdr3_group_per_sample,levels=rev(c('<0.01','[0.01,0.1)','[0.1,1)','[1,10)','>10')))

data@meta.data=meta

data_sub=subset(data,cells=rownames(meta)[which(!is.na(meta$vdj_barcode))])



meta_sub=meta[which(!is.na(meta$vdj_barcode)),]
c("#A6CEE3","royalblue","gold","brown")
meta_sub=cbind(data_sub@reductions$umap@cell.embeddings,meta_sub)

features=c('<0.01','[0.01,0.1)','[0.1,1)','[1,10)','>10')
plot_umap<-function(meta_sub,timepoint,res){
  tp=c('Baseline','On Induction')
  p=ggplot()+geom_point(meta_sub[which(meta_sub$clone_size_cdr3_group_per_sample==features[2]&meta_sub$Timepoint==timepoint&meta_sub$Response==res),],mapping=aes(x=UMAP_1,y=UMAP_2),size=0.1,color='#A6CEE3')+
    geom_point(meta_sub[which(meta_sub$clone_size_cdr3_group_per_sample==features[3]&meta_sub$Timepoint==timepoint&meta_sub$Response==res),],mapping=aes(x=UMAP_1,y=UMAP_2),size=0.1,color='royalblue')+
    geom_point(meta_sub[which(meta_sub$clone_size_cdr3_group_per_sample==features[4]&meta_sub$Timepoint==timepoint&meta_sub$Response==res),],mapping=aes(x=UMAP_1,y=UMAP_2),size=0.1,color='gold')+
    geom_point(meta_sub[which(meta_sub$clone_size_cdr3_group_per_sample==features[5]&meta_sub$Timepoint==timepoint&meta_sub$Response==res),],mapping=aes(x=UMAP_1,y=UMAP_2),size=0.1,color='brown')+
    theme_classic(base_size = 15)+theme(axis.title = element_blank(),axis.text = element_blank(),axis.line = element_blank(),axis.ticks =element_blank())+ggtitle(paste0(tp[timepoint],'-',res))
  return(p)
}

p1=plot_umap(meta_sub,1,'MHR')
p2=plot_umap(meta_sub,2,'MHR')

p3=plot_umap(meta_sub,1,'nonMHR')
p4=plot_umap(meta_sub,2,'nonMHR')

(p1+p2)/(p3+p4)

