library(philentropy)
library(ggplot2)
library(ggpubr)



set.seed(888)
######################################################Figure 6a  Transition index between MHR and nMHR###########################
#This code needs dataset file including patient id, Timepoint, cell subset annotation, 
#clonotype id or cdr3a+cdr3b sequence, whether patient have paired timepoints for each CD8+ T cell.
#All of the information and data can be download through XXXXX.
meta=readRDS('./cd8_meta_for_transition_index.rds')
#

meta$clonotype <- paste0(meta$PatientID,'_',meta$`CDR3A+CDR3B`)

pts=unique(meta$PatientID[which(meta$Paired!=0)])
df=data.frame(patient=pts)
for(i in 1:nrow(df)){
  kl=0
  temp=meta[which(meta$PatientID==pts[i]),]
  clones=unique(temp$clonotype)
  for(j in 1:length(clones)){
    temp.1=temp[which(temp$clonotype==clones[j]),]
    temp.2=table(temp.1$Timepoint)
   length(temp.2)
    if(length(temp.2)>1){
      temp.3=table(temp.1$Timepoint,temp.1$CellSubet)
      temp_pre=temp.3[1,]
      temp_post=temp.3[2,]
      kl.temp=KL(rbind(temp_pre,temp_post),est.prob = 'empirical')
      p.temp=temp.2['1']/length(which(temp$Timepoint==1))
      kl=kl+kl.temp*p.temp
    }
  }
    df[i,'Transition']=kl
    df[i,'Response']=unique(temp$Response)
}
df$Response[which(df$Response=='ER')]='MHR'
df$Response[which(df$Response=='R')]='nMHR'

ggplot(df[which(!is.na(df$Response)),],aes(x=Response,y=Transition,fill=Response))+geom_boxplot()+
  stat_compare_means()+geom_point()+theme_bw()+scale_fill_manual(values=c('royalblue','lightblue'))
