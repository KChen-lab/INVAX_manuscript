
library(stringr)
library(ggplot2)
library(rstatix)
library(ggpubr)

####{r clonotype statistics - indexes using whole vdj data}

get_tcr_stats<-function(patient_data,timepoint){
  patient_data_pre=patient_data[which(patient_data$Timepoint==timepoint),]
  temp=table(patient_data_pre$pseudo_clonotype,patient_data_pre$proportion_cdr3.aa_local  )
  temp=as.data.frame(temp)
  temp=temp[which(temp$Freq!=0),]
  temp$Var2=as.numeric(as.character(temp$Var2))
  
  if(nrow(temp)>0){
    shannon=-sum(temp$Var2 * log2(temp$Var2))
    evenness=shannon/log2(nrow(temp))  
    simpson=sum((temp$Var2)^2)
    simpson_div=1/sum((temp$Var2)^2)
    gini_simpson=1-simpson
    simpson_clonality=sqrt(simpson)
    return(c(shannon,evenness,simpson,simpson_div,gini_simpson,simpson_clonality))
  }
  else{
    return(c(NA,NA,NA,NA,NA,NA))
  }
  
}

wholevdj.dat=readRDS('./VDJ_data.rds')

response=read_xlsx('../DataShare/response.xlsx')

response=as.data.frame(response)
rownames(response)=response$Patient


df=data.frame(pts=str_replace(unique(wholevdj.dat$INVAX.Study.ID),'INVAX0',''),invax=unique(wholevdj.dat$INVAX.Study.ID))

df$response=as.character(response[match(as.numeric(df$pts),response$Patient),'Response'])

for(i in 1:nrow(df)){
  df[i,c('shannon','evenness','simpson','simpson_div','gini','simpson_clonolity')]=get_tcr_stats(wholevdj.dat[which(wholevdj.dat$INVAX.Study.ID==df$invax[i]),],1)
}

colnames(df)
df1=data.frame(invax=rep(df$invax,6),response=rep(df$response,6),value=c(df$shannon,df$evenness,df$simpson,df$simpson_div,df$gini,df$simpson_clonolity),
               index=c(rep('Shannon Entropy',nrow(df)),rep('Evenness',nrow(df)),rep('Simpson Index',nrow(df)),
                       rep('Simpson Diversity Index',nrow(df)), rep('Gini Simpson Index',nrow(df)), rep('Simpson clonality Index',nrow(df))))
df1=df1[which(df1$invax!='INVAX026'&df1$invax!='INVAX033'&!is.na(df1$response)),]
df1$response=as.character(df1$response)
df1$response[which(df1$response=='ER')]='MHR'
df1$response[which(df1$response=='R')]='nMHR'

p1=ggplot(df1[which(df1$index=='Shannon Entropy'),],aes(x=response,y=value,fill=response))+geom_boxplot() +
  geom_point()+stat_compare_means()+theme_classic(base_size = 15)+
  scale_fill_manual(values=c("darkblue","skyblue"))+ylab('Shannon Entropy')+
  xlab('')+theme(legend.position = 'none')
p2=ggplot(df1[which(df1$index=='Simpson Diversity Index'),],aes(x=response,y=value,fill=response))+geom_boxplot() +
  geom_point()+stat_compare_means()+theme_classic(base_size = 15)+
  scale_fill_manual(values=c("darkblue","skyblue"))+ylab('Simpson Diversity Index')+
  xlab('')+theme(legend.position = 'none')
p1+p2

df1=df1[which(df1$index=='Shannon Entropy'|df1$index=='Simpson Diversity Index'),]

stat.test=rbind(wilcox_test(value~response,data=df1[which(df1$index=='Shannon Entropy'),]),
                wilcox_test(value~response,data=df1[which(df1$index=='Simpson Diversity Index'),]))

