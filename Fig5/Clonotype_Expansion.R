library(rstatix)

#################################The rds file only contained annotation for expansion, you can try to run the calculation again using this code ##################
wholevdj.dat=readRDS('~/Documents/ChenLab/Projects/HeadNeck/final.figures/Manuscript/Data and code/SeuratObjects for reviewers/VDJ_data.rds')

raw_tcr_count=table(sapply(strsplit(wholevdj.dat$barcode,'_'),'[',1))
tt=table(wholevdj.dat$INVAX.Study.ID,wholevdj.dat$Timepoint)
pts=names(which(tt[,1]*tt[,2]!=0))#only consider paired samples


for (i in 1:length(pts)){
  temp=wholevdj.dat[which(wholevdj.dat$INVAX.Study.ID==pts[i]),]
  sample_sub=table(temp$Sample,temp$Timepoint)
  pre_sample=rownames(sample_sub)[which(sample_sub[,'1']!=0)]
  after_sample=rownames(sample_sub)[which(sample_sub[,'2']!=0)]

  temp_clone=unique(temp$pseudo_clonotype[!is.na(temp$pseudo_clonotype)])
  temp_count=table(temp$Timepoint)
  for(j in 1:length(temp_clone)){
    temp_temp=temp[which(temp$pseudo_clonotype==temp_clone[j]),]
    temp_count_temp=table(temp_temp$Timepoint)
    if(is.na(temp_count_temp['1'])){
      pre_freq=0
    }
    if(!is.na(temp_count_temp['1'])){
      pre_freq=unname(temp_count_temp['1']/temp_count['1'])
    }
    if(is.na(temp_count_temp['2'])){
      after_freq=0
    }
    if(!is.na(temp_count_temp['2'])){
      after_freq=unname(temp_count_temp['2']/temp_count['2'])
    }
    if(pre_freq==0){
      conti_table=rbind(c(0,temp_count['1']-0),
                        c(temp_count_temp['2'],temp_count['2']-temp_count_temp['2']))
    }
    if(after_freq==0){
      conti_table=rbind(c(temp_count_temp['1'],temp_count['1']-temp_count_temp['1']),
                        c(0,temp_count['2']-0))
    }
    if(pre_freq!=0&after_freq!=0){
      conti_table=rbind(c(temp_count_temp['1'],temp_count['1']-temp_count_temp['1']),
                        c(temp_count_temp['2'],temp_count['2']-temp_count_temp['2']))
    }
    conti_table=conti_table[c(2,1),]
    fishtest=fisher.test(conti_table)
    OR=unname(unlist(fishtest[3]))
    pval=unname(unlist(fishtest[1]))
    CI_low=fishtest$conf.int[1]
    CI_up=fishtest$conf.int[2]
    ind=which(wholevdj.dat$INVAX.Study.ID==pts[i]&wholevdj.dat$pseudo_clonotype==temp_clone[j])
    wholevdj.dat[ind,'pre_freq']=pre_freq
    wholevdj.dat[ind,'after_freq']=after_freq
    wholevdj.dat[ind,'Expansion_OR']=OR
    wholevdj.dat[ind,'Expansion_OR_pval']=pval
    wholevdj.dat[ind,'Expansion_OR_CI_low']=CI_low
    wholevdj.dat[ind,'Expansion_OR_CI_up']=CI_up
  }
}


pts=unique(wholevdj.dat$INVAX.Study.ID)
for(i in 1:length(pts)){
  patient_data=wholevdj.dat[which(wholevdj.dat$INVAX.Study.ID==pts[i]),]
  tt=table(patient_data$pseudo_clonotype,patient_data$Expansion_OR_pval)
  tt=as.data.frame(tt)
  tt=tt[which(tt$Freq!=0&tt$Freq!=1),]
  for( j in 1:nrow(tt)){
    tt[j,'pre']=length(which(patient_data$pseudo_clonotype==tt$Var1[j]&patient_data$Timepoint==1))
    tt[j,'post']=length(which(patient_data$pseudo_clonotype==tt$Var1[j]&patient_data$Timepoint==2))
  }
  temp=which(tt$pre==1&tt$post==1)
  tt=tt[-temp,] 
  tt$padj=p.adjust(as.numeric(as.character(tt$Var2)),method = 'fdr')
  for(j in 1:nrow(tt)){
    ind1=which(wholevdj.dat$INVAX.Study.ID==pts[i]&wholevdj.dat$pseudo_clonotype==as.character(tt$Var1[j]))
    wholevdj.dat[ind1,'Expansion_OR_pval_adj']=tt$padj[j]
  }
}

wholevdj.dat$Expansion_annot=NA
wholevdj.dat$Expansion_annot[which(!is.na(wholevdj.dat$Expansion_OR))]='None'
wholevdj.dat$Expansion_annot[which(wholevdj.dat$Expansion_OR<1&wholevdj.dat$Expansion_OR!=0&wholevdj.dat$Expansion_OR_pval_adj<0.05)]='Contracted_Sig'
wholevdj.dat$Expansion_annot[which(wholevdj.dat$Expansion_OR>1&wholevdj.dat$Expansion_OR_pval_adj<0.05)]='Expanded_Sig'
wholevdj.dat$Expansion_annot[which(wholevdj.dat$Expansion_OR==Inf&wholevdj.dat$Expansion_OR_pval_adj<0.05)]='Novel_Sig' 






colnames(wholevdj.dat)
table(wholevdj$Expansion_adj.v2)
table(wholevdj.dat$Expansion_annot)

meta_t=rbind(meta_cd4[,intersect(colnames(meta_cd4),colnames(meta_cd8))],meta_cd8[,intersect(colnames(meta_cd4),colnames(meta_cd8))])


setdiff(rownames(meta_t)[which(!is.na(meta_t$WholeVdj_barcode))],wholevdj$barcode)
ind1=match(intersect(wholevdj.dat$barcode,rownames(meta_t)),wholevdj.dat$barcode)
ind2=match(intersect(wholevdj.dat$barcode,rownames(meta_t)),rownames(meta_t))
identical(wholevdj.dat$barcode[ind1],rownames(meta_t)[ind2])
meta_t$WholeVdj_Expansion_adj.v2=NA
meta_t[ind2,'WholeVdj_Expansion_adj.v2']=wholevdj.dat[ind1,'Expansion_annot']


##from cd8 metadata, fraction against cd4+cd8;use this one
meta=meta_cd8
meta_t=rbind(meta_cd4[,intersect(colnames(meta_cd4),colnames(meta_cd8))],meta_cd8[,intersect(colnames(meta_cd4),colnames(meta_cd8))])
ind1=match(intersect(wholevdj.dat$barcode,rownames(meta_t)),wholevdj.dat$barcode)
ind2=match(intersect(wholevdj.dat$barcode,rownames(meta_t)),rownames(meta_t))
identical(wholevdj.dat$barcode[ind1],rownames(meta_t)[ind2])
meta_t$WholeVdj_Expansion_adj.v2=NA
meta_t[ind2,'WholeVdj_Expansion_adj.v2']=wholevdj.dat[ind1,'Expansion_annot']
ind1=match(intersect(wholevdj.dat$barcode,rownames(meta)),wholevdj.dat$barcode)
ind2=match(intersect(wholevdj.dat$barcode,rownames(meta)),rownames(meta))
identical(wholevdj.dat$barcode[ind1],rownames(meta)[ind2])
meta$WholeVdj_Expansion_adj.v2=NA
meta[ind2,'WholeVdj_Expansion_adj.v2']=wholevdj.dat[ind1,'Expansion_annot']

pts=unique(meta$meta.INVAX.Study.ID)
df=data.frame(invax=pts)
for(i in 1:nrow(df)){
  patient_data=meta[which(meta$meta.INVAX.Study.ID==df$invax[i]),]
  patient_data_all=meta_t[which(meta_t$meta.INVAX.Study.ID==df$invax[i]),]
  df[i,'Response']=unique(meta$PATH.Response.new[which(meta$meta.INVAX.Study.ID==df$invax[i])])
  df[i,'expanded_count']=length(unique(patient_data$WholeVdj_cdr3s[which(patient_data$WholeVdj_Expansion_adj.v2=='Expanded_Sig')]))
  df[i,'contracted_count']=length(unique(patient_data$WholeVdj_cdr3s[which(patient_data$WholeVdj_Expansion_adj.v2=='Contracted_Sig')]))
  df[i,'novel_count']=length(unique(patient_data$WholeVdj_cdr3s[which(patient_data$WholeVdj_Expansion_adj.v2=='Novel_Sig')]))
  df[i,'expanded_fraction']=length(unique(patient_data$WholeVdj_cdr3s[which(patient_data$WholeVdj_Expansion_adj.v2=='Expanded_Sig')]))/length(unique(patient_data_all$WholeVdj_cdr3s[which(patient_data_all$Timepoint==1)]))*100
  df[i,'contracted_fraction']=length(unique(patient_data$WholeVdj_cdr3s[which(patient_data$WholeVdj_Expansion_adj.v2=='Contracted_Sig')]))/length(unique(patient_data_all$WholeVdj_cdr3s[which(patient_data_all$Timepoint==1)]))*100
  df[i,'novel_fraction']=length(unique(patient_data$WholeVdj_cdr3s[which(patient_data$WholeVdj_Expansion_adj.v2=='Novel_Sig')]))/length(unique(patient_data_all$WholeVdj_cdr3s[which(patient_data_all$Timepoint==2)]))*100
}
df$Response[which(df$Response=='ER')]='MHR'
df$Response[which(df$Response=='R')]='nonMHR'
df=df[which(df$invax!='INVAX026'&df$invax!='INVAX033'),]
p1=ggplot(df[which(!is.na(df$Response)),],aes(x=Response,y=expanded_fraction,fill=Response))+geom_boxplot()+geom_point()+
  stat_compare_means()+scale_fill_manual(values=c('royalblue','lightblue'))+theme_classic(base_size = 15)+
  theme(legend.position = 'none',axis.text=element_text(color='black'))+ylab('% Fraction of expanded clones')+xlab('')
p2=ggplot(df[which(!is.na(df$Response)),],aes(x=Response,y=contracted_fraction,fill=Response))+geom_boxplot()+geom_point()+
  stat_compare_means()+scale_fill_manual(values=c('royalblue','lightblue'))+theme_classic(base_size = 15)+
  theme(legend.position = 'none',axis.text=element_text(color='black'))+ylab('% Fraction of contracted clones')+xlab('')
p3=ggplot(df[which(!is.na(df$Response)),],aes(x=Response,y=novel_fraction,fill=Response))+geom_boxplot()+geom_point()+
  stat_compare_means()+scale_fill_manual(values=c('royalblue','lightblue'))+theme_classic(base_size = 15)+
  theme(legend.position = 'none',axis.text=element_text(color='black'))+ylab('% Fraction of novel clones')+xlab('')
p1+p2+p3
p1+p2

ggsave('/Volumes/HTS-Intel-4/Research/HeadandNeck/hn.paper/final.figures/Figure5_final/5. Figure5CD_expanded_contracted_clones.pdf',
       width = 7,height = 4)
writexl::write_xlsx(df,'/Volumes/HTS-Intel-4/Research/HeadandNeck/hn.paper/final.figures/Figure5_final/SourceData/SourceData_5. Figure5CD_expanded_contracted_clones.xlsx')



