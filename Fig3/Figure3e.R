###############################################Figure 3e DEG of CD16+ and CD16- NK cells######################################

library(Seurat)
library(viridis)
library(ggrepel)
library(ggplot2)
library(ggpubr)
plot_volcano<-function(data,n1,n2){
  de=data
  fc_th=0.5
  de$DEG='NO'
  de$DEG[which(de$avg_log2FC>fc_th&de$p_val_adj<0.05)]='UP'
  de$DEG[which(de$avg_log2FC<(-fc_th)&de$p_val_adj<0.05)]='DOWN'
  if(length(which(de$DEG=='UP'))<20){
    de$DEG[order(de$avg_log2FC,decreasing = TRUE)[1:20]]='UP'
  }
  if(length(which(de$DEG=='DOWN'))<20){
    de$DEG[order(de$avg_log2FC,decreasing = FALSE)[1:20]]='DOWN'
  }
  de$delabel=NA
  de$delabel[which(de$DEG!='NO')]=rownames(de)[which(de$DEG!='NO')]
  de1=de[which(de$DEG=='UP'),]
  de2=de[which(de$DEG=='DOWN'),]
  de3=de[which(de$DEG=='NO'),]
  
  colors=viridis(3,option = 'A')
  
  p=ggplot(de,aes(x=avg_log2FC,y=-log10(p_val_adj),label=delabel))+
    geom_hline(yintercept = -log10(0.05),linetype="dotdash",color='brown',linewidth=1)+
    geom_vline(xintercept = -fc_th,linetype="dotdash",color='brown',linewidth=1)+
    geom_vline(xintercept = fc_th,linetype="dotdash",color='brown',linewidth=1)+
    geom_point(data=de3,aes(size=0.1),color='grey')+
    geom_point(data=de1,aes(x=avg_log2FC,y=-log10(p_val_adj),label=delabel,size=pct.2,color=pct.1))+
    geom_point(data=de2,aes(x=avg_log2FC,y=-log10(p_val_adj),label=delabel,size=pct.2,color=pct.1))+
    scale_color_viridis(option = 'D',paste0('Percent\nin ',n1))+
    scale_size_continuous(name=paste0('Percent\nin ',n2))+
    annotate(geom="text", x=-0.75, y=10, label=n2,color="blue",size=5)+
    annotate(geom="text", x=0.75, y=10, label=n1,color="red",size=5)+
    theme_bw(base_size = 15)+
    geom_text_repel(max.overlaps = 40)+
    ylab('-Log10(adjusted p value)')+xlab('Log2 fold change')
  return (p)
}



#################DEG can be recalculated using raw NK data downloaded from XXXXXX
nk=readRDS('./nk.rds') #This code need NK cell suerat object
####
meta_nk=nk@meta.data
nk_cd16=subset(nk,subset=FCGR3A>0)
meta_nk$FCGR3A_exp='None'
meta_nk[colnames(nk_cd16),'FCGR3A_exp']='CD16'
nk@meta.data=meta_nk
genes_de=grep(rownames(nk),pattern ="TRAV|TRAJ|TRAC|TRBV|TRBJ|TRBC|TRBD|TRDV|TRDJ|TRDC|TRDD|TRGV|TRGJ|TRGC|IGHV|IGLV|IGKV|IGHJ|IGLJ|IGKJ|IGHC|IGLC|IGKC|IGHE|IGHM|IGHA|IGHG|IGHD|IGHG|IGHA|IGHM|HPV16|JCHAIN|MALAT1|^RP([0-9]+-|[LS])|^HB[^(P)]",value = T )
nk=subset(nk,features = setdiff(rownames(nk),genes_de))
Idents(nk)='FCGR3A_exp'
cd16marker=FindMarkers(nk,ident.1 = 'CD16',ident.2='None',logfc.threshold = 0,min.pct = 0)
cd16marker$gene=rownames(cd16marker)
plot_volcano(cd16marker,'positive','negative')


#################DEG can be direcetly plotted using the supplementary table 9
#get the DEG list
cd16marker.deg=read_xlsx('../Supplementary Tables 6-13_v3.xlsx',
                sheet='Supplementary Table 9',range = 'A2:F99')
cd16marker.deg=as.data.frame(cd16marker.deg)
colnames(cd16marker.deg)=c('p_val','avg_log2FC','pct.1','pct.2','p_val_adj','gene')
rownames(cd16marker.deg)=cd16marker.deg$gene
plot_volcano(cd16marker.deg,'positive','negative')


