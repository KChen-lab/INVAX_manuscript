
library(ggalluvial)
library(ggplot2)

#########################Figure 5 halluvial plot
colnames(meta)

get_data <- function(metadata,expansion_type,response){
  meta=metadata[which(metadata$Expansion_annot==expansion_type&metadata$Response==response),]
  meta$group=meta$pseudoclonotype
  
  meta$Timepoint[which(meta$Timepoint==1)]='Pre'
  meta$Timepoint[which(meta$Timepoint==2)]='Post'
  
  states=rev(levels(meta$CellSubset))
  #states=states[c(1:11,13,14)]
  
  meta_temp=meta[which(!is.na(meta$vdj_barcode)),]
  #meta_temp=meta[which(!is.na(meta$WholeVdj_barcode)&meta$WholeVdj_frequency_cdr3.aa_local!=1),]
  
  #meta_temp=meta[which(!is.na(meta$WholeVdj_barcode)&meta$PATH.Response.new=='ER'),]
  count=0
  sk=data.frame()
  sk_heatmap=data.frame()
  for(i in 1:length(states)){
    for(j in 1:length(states)){
      count=count+1
      sk[count,'State1']=states[i]
      sk[count,'State2']=states[j]
      sk[count,'Freq_clone']=length(unique(intersect(meta_temp$group[which(meta_temp$CellSubset==states[i]&meta_temp$Timepoint=='Pre')],
                                                     meta_temp$group[which(meta_temp$CellSubset==states[j]&meta_temp$Timepoint=='Post')])))
      temp=unique(intersect(meta_temp$group[which(meta_temp$CellSubset==states[i]&meta_temp$Timepoint=='Pre')],
                            meta_temp$group[which(meta_temp$CellSubset==states[j]&meta_temp$Timepoint=='Post')]))
      sk[count,'Freq_clone']=length(temp)
      sk[count,'Freq_cell_pre']= length(which(!is.na(match(meta_temp$group[which(meta_temp$CellSubset==states[i]&meta_temp$Timepoint=='Pre')],temp))))
      sk[count,'Freq_cell_post']= length(which(!is.na(match(meta_temp$group[which(meta_temp$CellSubset==states[j]&meta_temp$Timepoint=='Post')],temp))))
    }
  }

  return(sk)
}




#This code needs:
#scTCR-seq dataset includes cell barcode pseudoclonotype,clonotype (defined by CDR3A+CDR3B) and 
#whether clonotypes has expanded or not; the expansion test is done in Clonotyope_Expansion.R 
#cd8 meta data includes patient id, response,timepoint

#data=readRDS('~/Documents/ChenLab/Projects/HeadNeck/final.figures/Manuscript/Data and code/SeuratObjects for reviewers/CD8.rds')
wholevdj.dat=readRDS('./VDJ_data.rds')
meta=readRDS('./CD8_metadata.rds')


#map scTCR data to cd8 data
ind1=match(intersect(wholevdj.dat$barcode,rownames(meta)),wholevdj.dat$barcode)
ind2=match(intersect(wholevdj.dat$barcode,rownames(meta)),rownames(meta))
identical(wholevdj.dat$barcode[ind1],rownames(meta)[ind2])
meta$Expansion_annot=NA
meta[ind2,'Expansion_annot']=wholevdj.dat[ind1,'Expansion_annot']
meta[ind2,'pseudoclonotype']=wholevdj.dat$pseudo_clonotype[ind1]
meta[ind2,'vdj_barcode']=wholevdj.dat$barcode[ind1]


type='Expanded_Sig'


colors=c("#A6CEE3", "#4F96C4","royalblue" ,"#CC3311", "#EE3377", "#FB9A99", "#A7D78D","#69BB54","#5D9E43" ,"#117733" , "#DDCC77", "#D9A398", "#CAB2D6" ,"#9A77B8", "#6A3D9A"  )  
states=levels(meta$CellSubset)
names(colors)=states

cells=states[c(2,3,4,6,7,8,10,11)]

sk=get_data(meta,type,'MHR')
ind=intersect(which(!is.na(match(sk$State1,cells))),which(!is.na(match(sk$State2,cells))))
sk=sk[ind,]
sk3.1=data.frame(State_trans=c(paste0(sk$State1,'_',sk$State2),paste0(sk$State1,'_',sk$State2)),
                 Freq=c(sk$Freq_cell_pre,sk$Freq_cell_post),timepoint=c(rep('Baseline',nrow(sk)),rep('OnTreatment',nrow(sk))),
                 State=c(sk$State1,sk$State2))
sk3.1$State=factor(sk3.1$State,levels=states[c(2,3,4,6,7,8,10,11)])
sk3.1$Response='MHR'

sk=get_data(meta,type,'nonMHR')
ind=intersect(which(!is.na(match(sk$State1,cells))),which(!is.na(match(sk$State2,cells))))
sk=sk[ind,]
sk3.2=data.frame(State_trans=c(paste0(sk$State1,'_',sk$State2),paste0(sk$State1,'_',sk$State2)),
                 Freq=c(sk$Freq_cell_pre,sk$Freq_cell_post),timepoint=c(rep('Baseline',nrow(sk)),rep('OnTreatment',nrow(sk))),
                 State=c(sk$State1,sk$State2))
sk3.2$State=factor(sk3.2$State,levels=states[c(2,3,4,6,7,8,10,11)])
sk3.2$Response='nonMHR'


sk3=rbind(sk3.1,sk3.2)


ggplot(sk3, aes(x = timepoint, stratum = State, alluvium = State_trans,y=Freq,
                fill = State, label = State)) +facet_wrap(.~Response)+
  geom_flow(stat = "alluvium", lode.guidance = "frontback") +
  geom_stratum()+theme_classic()+scale_fill_manual(values = colors[c(2,3,4,6,7,8,10,11)])+ylab('No of shared cell pairs')+xlab('')+theme_classic(base_size = 15)+
  ggtitle(paste0(type))


#########This code can be used to generate extended figure 8e