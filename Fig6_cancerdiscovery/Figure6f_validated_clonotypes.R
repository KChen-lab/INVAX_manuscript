
###################################################Figure 6f - Validated clonotypes###################################


library(dplyr)
library(Seurat)
library(ggplot2)


mycolors=c("#A6CEE3", "#4F96C4","royalblue" ,"#CC3311", "#EE3377", "#FB9A99", "#A7D78D","#69BB54","#5D9E43" ,"#117733" , "#DDCC77", "#D9A398", "#CAB2D6" ,"#9A77B8", "#6A3D9A"  )   
                                  

#this code need meta data for CD8 T cells, containing UMAP coordinates patient, cell subset, clonotype id for each cell
#all data can be downloaded from XXXXX
meta.1=readRDS('./cd8_meta_for_validated_clonotypes.rds')

#INVAX011	clonotype378 at baseline,clonotype1 on induction
#INVAX010	clonotype4 at baseline, clonotype1 on induction
#INVAX021	clonotype71 at baseline, clonotype8 on induction

#three validated clonotypes
ind=list(clone1=c(which(meta.1$Patient=='INVAX011'&meta.1$RawClonotypeID=='clonotype378'&meta.1$Timepoint==1),
                  which(meta.1$Patient=='INVAX011'&meta.1$RawClonotypeID=='clonotype1'&meta.1$Timepoint==2)),
         
         clone2=c(which(meta.1$Patient=='INVAX021'&meta.1$RawClonotypeID=='clonotype71'&meta.1$Timepoint==1),
                  which(meta.1$Patient=='INVAX021'&meta.1$RawClonotypeID=='clonotype8'&meta.1$Timepoint==2)),
         clone3=c(which(meta.1$Patient=='INVAX010'&meta.1$RawClonotypeID=='clonotype4'&meta.1$Timepoint==1),
                  which(meta.1$Patient=='INVAX010'&meta.1$RawClonotypeID=='clonotype1'&meta.1$Timepoint==2))
         )


for(i in 1:3){
df=meta.1
ttt=as.character(df$CellSubset)
df$Highlight_TCR=0
df$Highlight_TCR[ind[[i]]]='Baseline'
df$Highlight_TCR[which(df$Highlight_TCR=='Baseline'&df$Timepoint==2)]=ttt[which(df$Highlight_TCR=='Baseline'&df$Timepoint==2)]
df$Highlight_TCR=factor(df$Highlight_TCR,levels=c(levels(df$CellSubset),0,'Baseline'))

df$Timepoint=factor(df$Timepoint)
#plot UMAP
ggplot()+ geom_point(df[which(df$Highlight_TCR==0),], mapping=aes(UMAP_1, UMAP_2), size=0.1,color='grey') + 
  # all other cells are in grey
  geom_point(df[which(str_detect(df$Highlight_TCR,'Baseline')),], mapping=aes(UMAP_1, UMAP_2,color=Highlight_TCR), size=2,color='black')+
  geom_point(df[which(df$Highlight_TCR!=0&!str_detect(df$Highlight_TCR,'Baseline')),], mapping=aes(UMAP_1, UMAP_2,color=Highlight_TCR), size=2)+ 
  scale_color_manual(values=c(mycolors[which(!is.na(match(levels(meta.1$CellSubset),unique(df$Highlight_TCR))))]))+
  #cells from baseline is black
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),legend.title = element_blank(),
        axis.ticks.y=element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
        axis.title = element_blank())+ggtitle(paste0('Clonotype ',i))
ggsave(paste0('./UMAP_validated_clonotype.',i,'.pdf'),height = 4,width = 6)

#plot barplot
temp=df[ind[[i]],]
tt=as.data.frame.matrix(table(temp$Timepoint,temp$CellSubset))
tt=as.data.frame(table(temp$Timepoint,temp$CellSubset))
tt$Var1=as.numeric(tt$Var1)
tt$Var1[which(tt$Var1==1)]='Baseline'
tt$Var1[which(tt$Var1==2)]='OnTreatment'
colnames(tt)=c('Timepoint','CD8Subset','Count')
ggplot(tt,aes(x=CD8Subset,y=Count,fill=Timepoint))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_classic()+scale_fill_manual(values = c("#DDCC77","#AA4499"))+theme(axis.text.x = element_text(angle = 90,size=10,hjust = 1,vjust = 0.5))
ggsave(paste0('./Barplot_validated_clonotype.',i,'.pdf'),height = 4,width = 6)

}


                                  
            
                                  