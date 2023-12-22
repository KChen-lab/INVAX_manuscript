
library(Seurat)
library(ggplot2)

##############################Figure 5E umap of expanded and contracted clonotypes###########




#This code needs:
#scTCR-seq dataset includes cell barcode pseudoclonotype,clonotype (defined by CDR3A+CDR3B) and 
#whether clonotypes has expanded or not; the expansion test is done in Clonotyope_Expansion.R 
#cd8 meta data includes patient id, response,timepoint
data=readRDS('./CD8.rds')
wholevdj.dat=readRDS('./VDJ_data.rds')
#{r clonotype size umap }
meta=data@meta.data


#map scTCR data to cd8 data
ind1=match(intersect(wholevdj.dat$barcode,rownames(meta)),wholevdj.dat$barcode)
ind2=match(intersect(wholevdj.dat$barcode,rownames(meta)),rownames(meta))
identical(wholevdj.dat$barcode[ind1],rownames(meta)[ind2])
meta$Expansion_annot=NA
meta[ind2,'Expansion_annot']=wholevdj.dat[ind1,'Expansion_annot']
meta[ind2,'pseudoclonotype']=wholevdj.dat$pseudo_clonotype[ind1]
meta[ind2,'vdj_barcode']=wholevdj.dat$barcode[ind1]



type='Expanded_Sig'
meta[,c('UMAP_1','UMAP_2')]=data[["umap"]]@cell.embeddings
meta_temp2=meta[which(meta$Expansion_annot!=type),]
meta_temp1=meta[which(meta$Expansion_annot==type),]

colors=c("#A6CEE3", "#4F96C4","royalblue" ,"#CC3311", "#EE3377", "#FB9A99", "#A7D78D","#69BB54","#5D9E43" ,"#117733" , "#DDCC77", "#D9A398", "#CAB2D6" ,"#9A77B8", "#6A3D9A"  )  

cols=colors[match(intersect(levels(meta_temp1$CellSubset),unique(meta_temp1$CellSubset[which(meta_temp1$Timepoint==1&meta_temp1$Response=='MHR')])),levels(meta_temp1$CellSubset))]
p1=ggplot()+ geom_point(meta_temp2, mapping=aes(UMAP_1, UMAP_2),alpha=0.25,color='grey')+
  geom_point(meta_temp1[which(meta_temp1$Timepoint==1&meta_temp1$Response=='MHR'),],mapping=aes(UMAP_1, UMAP_2,color=CellSubset), size=1)+scale_color_manual(values=cols)+theme_classic(base_size = 20)+
  theme(legend.position = 'none',axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())+ggtitle('Baseline:MHR')
cols=colors[match(intersect(levels(meta_temp1$CellSubset),unique(meta_temp1$CellSubset[which(meta_temp1$Timepoint==2&meta_temp1$Response=='MHR')])),levels(meta_temp1$CellSubset))]
p2=ggplot()+ geom_point(meta_temp2, mapping=aes(UMAP_1, UMAP_2),alpha=0.25,color='grey')+
  geom_point(meta_temp1[which(meta_temp1$Timepoint==2&meta_temp1$Response=='MHR'),],mapping=aes(UMAP_1, UMAP_2,color=CellSubset), size=1)+scale_color_manual(values=cols)+theme_classic(base_size = 20)+theme(legend.position = 'none',axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())+ggtitle('OnTreatment:MHR')
cols=colors[match(intersect(levels(meta_temp1$CellSubset),unique(meta_temp1$CellSubset[which(meta_temp1$Timepoint==1&meta_temp1$Response=='nonMHR')])),levels(meta_temp1$CellSubset))]
p3=ggplot()+ geom_point(meta_temp2, mapping=aes(UMAP_1, UMAP_2),alpha=0.25,color='grey')+
  geom_point(meta_temp1[which(meta_temp1$Timepoint==1&meta_temp1$Response=='nonMHR'),],mapping=aes(UMAP_1, UMAP_2,color=CellSubset), size=1)+scale_color_manual(values=cols)+theme_classic(base_size = 20)+theme(legend.position = 'none',axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())+ggtitle('Baseline:nonMHR')
cols=colors[match(intersect(levels(meta_temp1$CellSubset),unique(meta_temp1$CellSubset[which(meta_temp1$Timepoint==2&meta_temp1$Response=='nonMHR')])),levels(meta_temp1$CellSubset))]
p4=ggplot()+ geom_point(meta_temp2, mapping=aes(UMAP_1, UMAP_2),alpha=0.25,color='grey')+
  geom_point(meta_temp1[which(meta_temp1$Timepoint==2&meta_temp1$Response=='nonMHR'),],mapping=aes(UMAP_1, UMAP_2,color=CellSubset), size=1)+scale_color_manual(values=cols)+theme_classic(base_size = 20)+theme(legend.position = 'none',axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())+ggtitle('OnTreatment:nonMHR')
(p1|p2)/(p3|p4)


##############################Figure 5f
df1=table(meta_temp1$CellSubset[which(meta_temp1$Response=='MHR')],meta_temp1$Timepoint[which(meta_temp1$Response=='MHR')])
df2=table(meta_temp1$CellSubset[which(meta_temp1$Response=='nonMHR')],meta_temp1$Timepoint[which(meta_temp1$Response=='nonMHR')])
df1=as.data.frame(df1)
df2=as.data.frame(df2)
df1$Response='MHR'
df2$Response='nonMHR'
df=rbind(df1,df2)
ggplot(df,aes(x=Var1,y=Freq,fill=Var2))+
  geom_bar(stat='identity',position = position_dodge())+facet_wrap(.~Response,ncol=1)+theme_classic(base_size = 15)+
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5),legend.position = 'none') +
  ylab('No. of cells')+
  scale_fill_manual(values=c('gold','purple'))+xlab('CD8 cluster')


############This code can be also used for extended figure 8b-c,8f-g by changing the type name