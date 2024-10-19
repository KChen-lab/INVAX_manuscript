
library(ggplot2)
library(rstatix)
library(ggpubr)


######{r 10 - 0.1% frequency clone number}

################################################This code needs cd4 and cd8 cell metadata############
wholevdj.dat=readRDS('./VDJ_data.rds')
meta=readRDS('./CD8_metadata.rds')
meta.cd4=readRDS('./CD4_metadata.rds')
meta_total=rbind(meta.cd4[,intersect(colnames(meta.cd4),colnames(meta))],meta[,intersect(colnames(meta.cd4),colnames(meta))])


temp=intersect(wholevdj.dat$barcode,rownames(meta_total))
ind1=match(temp,wholevdj.dat$barcode)
ind2=match(temp,rownames(meta_total))
identical(wholevdj.dat$barcode[ind1],rownames(meta_total)[ind2])
meta_total[ind2,'pseudoclonotype']=wholevdj.dat$pseudo_clonotype[ind1]
meta_total[ind2,'proportion_cdr3.aa_local']=wholevdj.dat$proportion_cdr3.aa_local[ind1]
meta_total[ind2,'vdj_barcode']=wholevdj.dat$barcode[ind1]

temp=intersect(wholevdj.dat$barcode,rownames(meta))
ind1=match(temp,wholevdj.dat$barcode)
ind2=match(temp,rownames(meta))
identical(wholevdj.dat$barcode[ind1],rownames(meta)[ind2])
meta[ind2,'pseudoclonotype']=wholevdj.dat$pseudo_clonotype[ind1]
meta[ind2,'proportion_cdr3.aa_local']=wholevdj.dat$proportion_cdr3.aa_local[ind1]
meta[ind2,'vdj_barcode']=wholevdj.dat$barcode[ind1]


meta=meta[which(meta$patient!='INVAX026'),]
meta_total=meta_total[which(meta_total$patient!='INVAX026'),]

freqs=seq(0.001,0.1,0.001)
t_plot=data.frame(freqs=freqs)
for (i in 1:length(freqs)){
  t_plot[i,'MHR']=length(unique(paste0(meta$pseudoclonotype,meta$patient)[which(meta$proportion_cdr3.aa_local>freqs[i]&meta$Timepoint==1&meta$Response=='MHR')]))
  t_plot[i,'nonMHR']=length(unique(paste0(meta$pseudoclonotype,meta$patient)[which(meta$proportion_cdr3.aa_local>freqs[i]&meta$Timepoint==1&meta$Response=='nonMHR')]))
  
}
t_plot$total_ER=length(unique(paste0(meta_total$pseudoclonotype,meta_total$patient)[which(meta_total$Timepoint==1&meta_total$Response=='MHR')]))
t_plot$total_R=length(unique(paste0(meta_total$pseudoclonotype,meta_total$patient)[which(meta_total$Timepoint==1&meta_total$Response=='nonMHR')]))

t_plot_1=data.frame(freqs=c(t_plot$freqs,t_plot$freqs),Response=c(rep('MHR',nrow(t_plot)),rep('nonMHR',nrow(t_plot))),
                    clone=c(t_plot$MHR,t_plot$nonMHR),clone_norm=c(t_plot$MHR/t_plot$total_ER,t_plot$nonMHR/t_plot$total_R))

t_plot_1$clone_norm=t_plot_1$clone_norm*100
ggplot(data=t_plot_1, aes(x=freqs, y=clone_norm, group=Response,color=Response)) +geom_vline(xintercept = 0.01,size=1,linetype=2)+
  geom_line(linetype = "dashed")+
  geom_point()+scale_color_manual(values=c( "royalblue","lightblue"))+scale_y_continuous(trans = 'log10')+
  theme_classic(base_size = 15)+ylab('Fration of CD8 clones (%)')+xlab('Clonotype frequency thresholds')+
  theme(axis.text=element_text(color='black'))



