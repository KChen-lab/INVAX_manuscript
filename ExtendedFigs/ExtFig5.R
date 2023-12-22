
###############################################Extended Figure5 NK cell clusters ####################################

library(Seurat)
library(ggplot2)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(rstatix)
nk=readRDS('./NK.rds')


##Extended Figure 5a
color_palette <- c(
  "#FF0000", # Fire Engine Red
  "#0000FF", # Electric Blue
  "#39FF14", # Neon Green
  "#FFD700", # Sunflower Yellow
  "#800080", # Deep Purple
  "#FFA500", # Tangerine Orange
  "#FF69B4", # Hot Pink
  "#008080", # Teal
  'salmon',
  "#696969"  # Charcoal Gray
)
DimPlot(nk, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6,cols = color_palette) +NoAxes()

##Extended Figure5b
m=c('IL16','IL10RA','IL10RB','CXCR6','IL2RB','XCL2','XCL1','GZMK','CXCR3','CCR1','CCR7','IL7R','IL18','IL18R1','IL4R','IL2RA',
    'CXCR4','IL21R','PRF1','GZMB','GZMH','CXCR2','IL2RG','TGFBR1','CX3CR1','TGFBR3','TGFBR2','GZMA','CCL5','IL32','CCL4L2','CCL4','CCL3')

Idents(nk)
exp=AverageExpression(nk,features = m)
pheatmap(t(exp$RNA),scale = 'column',cluster_cols = F,cluster_rows = F,name = 'avg.exp',color = colorRampPalette(c("navy", 'white',"red"))(50),border_color = 'white')


##########Extended Figure5d
meta_nk=nk@meta.data
scoresnk=meta_nk[,c('Cytotoxicity','Stimulatory','Inhibitory','Cyto2Inhib')]
df=data.frame(cluster=levels(meta_nk$CellSubset))
for (i in 1:nrow(df)){
  df[i,'Cytotoxicity']=mean(meta_nk$Cytotoxicity[which(meta_nk$CellSubset==df$cluster[i])])
  df[i,'Stimulatory']=mean(meta_nk$Stimulatory[which(meta_nk$CellSubset==df$cluster[i])])
  df[i,'Inhibitory']=mean(meta_nk$Inhibitory[which(meta_nk$CellSubset==df$cluster[i])])
  df[i,'Ratio']=mean(meta_nk$Cyto2Inhib[which(meta_nk$CellSubset==df$cluster[i])])
  
}
rownames(df)=df$cluster
df1=df[order(df$Cytotoxicity,decreasing = T),]
pdf('~/Documents/ChenLab/Projects/HeadNeck/final.figures/Figure 3/NK_cluster_scores_heatmap_202312.pdf',width = 3,height = 6)
pheatmap(df1[,c(3,2,4,5)],scale = 'column',cluster_rows = F,cluster_cols  = F,color = colorRampPalette(c("navy", "white", "red"))(50),main = 'Gene set score')
dev.off()
pheatmap(df1[,c(3,2,4,5)],scale = 'column',cluster_rows = F,cluster_cols  = F,color = colorRampPalette(c("navy", "white", "red"))(50),main = 'Gene set score')



##########Extended Figure5e
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
p1=VlnPlot(subset(nk,idents='NK.c5',subset=Timepoint==1&(Response=='MHR'|Response=='nonMHR')),features = c('Cytotoxicity'),group.by = 'Response',pt.size = 0,
           cols = c('royalblue','lightblue'))&stat_compare_means()&
  stat_summary(fun.data = data_summary,geom='pointrange', size = 0.5, colour = "yellow")&NoLegend()
p2=VlnPlot(subset(nk,idents='NK.c5',subset=Timepoint==1&(Response=='MHR'|Response=='nonMHR')),features = c('Stimulatory'),group.by = 'Response',pt.size = 0,
           cols = c('royalblue','lightblue'))&stat_compare_means()&
  stat_summary(fun.data = data_summary,geom='pointrange', size = 0.5, colour = "yellow")&NoLegend()
p3=VlnPlot(subset(nk,idents='NK.c5',subset=Timepoint==1&(Response=='MHR'|Response=='nonMHR')),features = c('FCGR3A'),group.by = 'Response',pt.size = 0,
           cols =c('royalblue','lightblue'))&stat_compare_means()&
  stat_summary(fun.data = data_summary,geom='pointrange', size = 0.5, colour = "yellow")&NoLegend()
p1=p1+xlab('')+ylab('Gene set score')
p2=p2+xlab('')+ylab('Gene set score')
p3=p3+xlab('')+ylab('Expression level')
p1|p2|p3

nk.sub=subset(nk,idents='NK.c5',subset=Timepoint==1&(Response=='MHR'|Response=='nonMHR'))
meta.nk.sub=nk.sub@meta.data
meta.nk.sub$FCGR3A=FetchData(nk.sub,vars = 'FCGR3A')[,1]
cyto.1=data_summary(nk.sub$Cytotoxicity[which(nk.sub$Response=='MHR')])
cyto.2=data_summary(nk.sub$Cytotoxicity[which(nk.sub$Response=='nonMHR')])
cyto=wilcox_test(Cytotoxicity~Response,data=meta.nk.sub)
stim.1=data_summary(nk.sub$Stimulatory[which(nk.sub$Response=='MHR')])
stim.2=data_summary(nk.sub$Stimulatory[which(nk.sub$Response=='nonMHR')])
stim=wilcox_test(Stimulatory~Response,data=meta.nk.sub)
FCGR3A.1=data_summary(meta.nk.sub$FCGR3A[which(meta.nk.sub$Response=='MHR')])
FCGR3A.2=data_summary(meta.nk.sub$FCGR3A[which(meta.nk.sub$Response=='nonMHR')])
FCGR3A=wilcox_test(FCGR3A~Response,data=meta.nk.sub)

df_stat=data.frame(value=c(meta.nk.sub$Cytotoxicity,meta.nk.sub$Stimulatory,meta.nk.sub$FCGR3A),
                   name=c(rep('Cytotoxicity',nrow(meta.nk.sub)),rep('Stimulatory',nrow(meta.nk.sub)),rep('FCGR3A',nrow(meta.nk.sub))),
                   response=c(meta.nk.sub$Response,meta.nk.sub$Response,meta.nk.sub$Response))
compare_means(value~response,data=df_stat,group.by = 'name')

stats=rbind(cyto.1,cyto.2,stim.1,stim.2,FCGR3A.1,FCGR3A.2)
rownames(stats)=c('Cytotoxcity.MHR','Cytotoxcity.nMHR',
                  'Stimulatory.MHR','Stimulatory.nMHR',
                  'FCGR3A.MHR','FCGR3A.nMHR')

stats.1=rbind(cyto,stim,FCGR3A)
