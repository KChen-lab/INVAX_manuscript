

library(reshape2)
library(scales)
library(cowplot)
library(ggplot2)
library(readxl)
library(dplyr)
library(rstatix)
library(ggpubr)

stat1=read_xlsx('~/Documents/ChenLab/Projects/HeadNeck/final.figures/Manuscript/1207/Supplementary Tables 6-13_v3.xlsx',
                sheet='Supplementary Table 7',range = 'A3:AE73')
response=read_xlsx('../DataShare/response.xlsx')
response=as.data.frame(response)

stat1$Response=as.character(response[match(stat1$`Patient ID`,response$Patient),'Response'])
stat1$Response=factor(stat1$Response,levels=c('MHR','nMHR'))

###########################Extended Figure 3c-d#################################
#comparing between response at baseline
stat1=stat1[which(!is.na(stat1$Response)),]
stat1

tt_MHR=data.frame(CellSubset = c(names(colSums(stat1[which(stat1$Timepoint=='Baseline'&stat1$Response=='MHR'&stat1$`Patient ID`!=33),][,3:31])),
                                 names(colSums(stat1[which(stat1$Timepoint=='OnInduction'&stat1$Response=='MHR'),][,3:31]))),
           Freq = c(unlist(colSums(stat1[which(stat1$Timepoint=='Baseline'&stat1$Response=='MHR'&stat1$`Patient ID`!=33),][,3:31])),
                    unlist(colSums(stat1[which(stat1$Timepoint=='OnInduction'&stat1$Response=='MHR'),][,3:31]))),
           Timepoint=c(rep(1,29),rep(2,29))
           )

tt_nonMHR=data.frame(CellSubset = c(names(colSums(stat1[which(stat1$Timepoint=='Baseline'&stat1$Response=='nMHR'),][,3:31])),
                                 names(colSums(stat1[which(stat1$Timepoint=='OnInduction'&stat1$Response=='nMHR'),][,3:31]))),
                  Freq = c(unlist(colSums(stat1[which(stat1$Timepoint=='Baseline'&stat1$Response=='nMHR'),][,3:31])),
                           unlist(colSums(stat1[which(stat1$Timepoint=='OnInduction'&stat1$Response=='nMHR'),][,3:31]))),
                  Timepoint=c(rep(1,29),rep(2,29))
)

cells=as.character(unique(tt_MHR$CellSubset))
tt_plot=data.frame(cells=as.character(unique(tt_MHR$CellSubset)))
for(i in 1:length(cells)){
  tt_plot[i,'pre_ER']=tt_MHR$Freq[which(tt_MHR$CellSubset==cells[i]&tt_MHR$Timepoint==1)]/sum(tt_MHR$Freq[which(tt_MHR$Timepoint==1)])
  tt_plot[i,'post_ER']=tt_MHR$Freq[which(tt_MHR$CellSubset==cells[i]&tt_MHR$Timepoint==2)]/sum(tt_MHR$Freq[which(tt_MHR$Timepoint==2)])
  tt_plot[i,'pre_NER']=tt_nonMHR$Freq[which(tt_nonMHR$CellSubset==cells[i]&tt_nonMHR$Timepoint==1)]/sum(tt_nonMHR$Freq[which(tt_nonMHR$Timepoint==1)])
  tt_plot[i,'post_NER']=tt_nonMHR$Freq[which(tt_nonMHR$CellSubset==cells[i]&tt_nonMHR$Timepoint==2)]/sum(tt_nonMHR$Freq[which(tt_nonMHR$Timepoint==2)])
  ind=which(tt_MHR$CellSubset==cells[i])
  temp=fisher.test(matrix(c(tt_MHR$Freq[ind[2]],
                            sum(tt_MHR$Freq[which(tt_MHR$Timepoint==2)])-tt_MHR$Freq[ind[2]],
                            tt_MHR$Freq[ind[1]],
                            sum(tt_MHR$Freq[which(tt_MHR$Timepoint==1)])-tt_MHR$Freq[ind[1]]),
                          nrow = 2))
  tt_plot[i,'pval_ER']=temp$p.value
  tt_plot[i,'odds_ER']=log2(temp$estimate)
  tt_plot[i,'CI.low_ER']=log2(temp$conf.int[1])
  tt_plot[i,'CI.high_ER']=log2(temp$conf.int[2])

  ind=which(tt_nonMHR$CellSubset==cells[i])
  temp=fisher.test(matrix(c(tt_nonMHR$Freq[ind[2]],
                            sum(tt_nonMHR$Freq[which(tt_nonMHR$Timepoint==2)])-tt_nonMHR$Freq[ind[2]],
                            tt_nonMHR$Freq[ind[1]],
                            sum(tt_nonMHR$Freq[which(tt_nonMHR$Timepoint==1)])-tt_nonMHR$Freq[ind[1]]),
                          nrow = 2))
  tt_plot[i,'pval_NER']=temp$p.value
  tt_plot[i,'odds_NER']=log2(temp$estimate)
  tt_plot[i,'CI.low_NER']=log2(temp$conf.int[1])
  tt_plot[i,'CI.high_NER']=log2(temp$conf.int[2])
}

tt_plot$pval_ER.adj=p.adjust(tt_plot$pval_ER,method = 'fdr')
tt_plot$pval_NER.adj=p.adjust(tt_plot$pval_NER,method = 'fdr')

tt_plot1=tt_plot[15:29,]
tt_plot2=data.frame(Abundance=c(tt_plot1$pre_ER,tt_plot1$pre_NER,tt_plot1$post_ER,tt_plot1$post_NER),
                    response=c(rep('MHR',nrow(tt_plot1)),rep('nonMHR',nrow(tt_plot1)),rep('MHR',nrow(tt_plot1)),rep('nonMHR',nrow(tt_plot1))),
                    cells=c(tt_plot1$cells,tt_plot1$cells,tt_plot1$cells,tt_plot1$cells),
                    pvals=c(as.character(symnum(c(tt_plot1$pval_ER.adj,tt_plot1$pval_NER.adj), corr = FALSE, na = FALSE, cutpoints = c(0,  0.001,  1), symbols = c( "*", " "))),rep(NA,2*nrow(tt_plot1))),
                    timepoint=c(rep('Baseline',2*nrow(tt_plot1)),rep('OnInduction',2*nrow(tt_plot1))),
                    odds=c(tt_plot1$odds_ER,tt_plot1$odds_NER,rep(NA,2*nrow(tt_plot1)))
)


tt_plot2$cells=factor(tt_plot2$cells,levels=rev(tt_plot$cells[15:29]))
tt_plot2$timepoint=factor(tt_plot2$timepoint,levels=c('OnInduction','Baseline'))
tt_plot2$Abundance=100*tt_plot2$Abundance

tt_plot2.1=tt_plot2[which(tt_plot2$response=='MHR'),]
or=order(tt_plot2.1$odds[1:15])
tt_plot2$cells=factor(tt_plot2.1$cells,levels=rev(levels(tt_plot2.1$cells))[or])

p1=ggplot(tt_plot2[which(tt_plot2$response=='MHR'),],aes(y=cells,x=Abundance,fill=timepoint,label=pvals))+geom_bar(stat = 'identity',position = position_dodge(width = 0.8),color='black')+facet_grid(.~response)+
  scale_fill_manual(values=c( "#AA4499","#DDCC77"),'Timepoint')+theme_classic(base_size = 15)+
  theme(axis.text.x = element_text(angle = 90,hjust = 0.95,vjust = 0.5,color='black'),axis.text.y = element_text(color='black'))+xlab('% in all T cells')+ylab('')+
  geom_text(size=8)

p2=ggplot(tt_plot2[1:15,],aes(y=cells,fill=odds,x=-0.5))+geom_tile(color='grey30')+scale_fill_gradientn(colours = c("navy", "white", "firebrick"),limits = c(min(tt_plot2$odds[1:30]), max(tt_plot2$odds[1:30])),
                                                                                                        values = rescale(c(min(tt_plot2$odds[1:30]),0, max(tt_plot2$odds[1:30]))),name='Log2\nodds\nratio')+
  theme_classic(base_size = 15 )+theme(axis.title.x = element_blank(),axis.text.x  = element_blank(),axis.ticks.x  = element_blank(),
                                       legend.position = 'left')
p3=ggplot(tt_plot2[which(tt_plot2$response=='nonMHR'),],aes(y=cells,x=Abundance,fill=timepoint,label=pvals))+geom_bar(stat = 'identity',position = position_dodge(width = 0.8),color='black')+facet_grid(.~response)+
  scale_fill_manual(values=c( "#AA4499","#DDCC77"),'Timepoint')+theme_classic(base_size = 15)+
  theme(axis.text.x = element_text(angle = 90,hjust = 0.95,vjust = 0.5,color='black'),axis.text.y = element_text(color='black'))+xlab('% in all T cells')+ylab('')+
  geom_text(size=8)
p4=ggplot(tt_plot2[16:30,],aes(y=cells,fill=odds,x=-0.5))+geom_tile(color='grey30')+scale_fill_gradientn(colours = c("navy", "white", "firebrick"),limits = c(min(tt_plot2$odds[1:30]), max(tt_plot2$odds[1:30])),
                                                                                                         values = rescale(c(min(tt_plot2$odds[1:30]),0, max(tt_plot2$odds[1:30]))),name='Log2\nodds\nratio')+
  theme_classic(base_size = 15 )+theme(axis.title.x = element_blank(),axis.text.x  = element_blank(),axis.ticks.x  = element_blank(),legend.position = 'left')
plot_grid(p2+ylab(''),
          p1+theme(axis.title.y = element_blank(),axis.text.y  = element_blank(),axis.ticks.y  = element_blank(),axis.line.y = element_blank(),legend.position = 'none')+ylab(''),
          p4+ylab(''),
          p3+theme(axis.title.y = element_blank(),axis.text.y  = element_blank(),axis.ticks.y  = element_blank(),axis.line.y = element_blank(),legend.position = 'none')+ylab(''),
          align = "h", ncol = 4, axis = "tb", rel_widths  = c(6, 10,6,10))




# This code can also be used to generate figure for extended figure 3e
# tt_plot1=tt_plot[1:14,]
# tt_plot2=data.frame(Abundance=c(tt_plot1$pre_ER,tt_plot1$pre_NER,tt_plot1$post_ER,tt_plot1$post_NER),
#                     response=c(rep('MHR',nrow(tt_plot1)),rep('nonMHR',nrow(tt_plot1)),rep('MHR',nrow(tt_plot1)),rep('nonMHR',nrow(tt_plot1))),
#                     cells=c(tt_plot1$cells,tt_plot1$cells,tt_plot1$cells,tt_plot1$cells),
#                     pvals=c(as.character(symnum(c(tt_plot1$pval_ER.adj,tt_plot1$pval_NER.adj), corr = FALSE, na = FALSE, cutpoints = c(0,  0.001,  1), symbols = c( "*", " "))),rep(NA,2*nrow(tt_plot1))),
#                     timepoint=c(rep('Baseline',2*nrow(tt_plot1)),rep('OnInduction',2*nrow(tt_plot1))),
#                     odds=c(tt_plot1$odds_ER,tt_plot1$odds_NER,rep(NA,2*nrow(tt_plot1)))
# )
# tt_plot2$cells=factor(tt_plot2$cells,levels=rev(tt_plot$cells[1:14]))
# tt_plot2$Abundance=100*tt_plot2$Abundance
# tt_plot2$timepoint=factor(tt_plot2$timepoint,levels=c('OnInduction','Baseline'))
# 
# 
# tt_plot2.1=tt_plot2[which(tt_plot2$response=='MHR'),]
# or=order(tt_plot2.1$odds[1:14])
# tt_plot2$cells=factor(tt_plot2.1$cells,levels=rev(levels(tt_plot2.1$cells))[or])
# p1=ggplot(tt_plot2[which(tt_plot2$response=='MHR'),],aes(y=cells,x=Abundance,fill=timepoint,label=pvals))+geom_bar(stat = 'identity',position = position_dodge(width = 0.8),color='black')+facet_grid(.~response)+
#   scale_fill_manual(values=c( "#AA4499","#DDCC77"),'Timepoint')+theme_classic(base_size = 15)+
#   theme(axis.text.x = element_text(angle = 90,hjust = 0.95,vjust = 0.5,color='black'),axis.text.y = element_text(color='black'))+xlab('% in all T cells')+ylab('')+
#   geom_text(size=8)
# p2=ggplot(tt_plot2[1:14,],aes(y=cells,fill=odds,x=-0.5))+geom_tile(color='grey30')+scale_fill_gradientn(colours = c("navy", "white", "firebrick"),limits = c(min(tt_plot2$odds[1:28]), max(tt_plot2$odds[1:28])),
#                                                                                                         values = rescale(c(min(tt_plot2$odds[1:28]),0, max(tt_plot2$odds[1:28]))),name='Log2\nodds\nratio')+
#   theme_classic(base_size = 15 )+theme(axis.title.x = element_blank(),axis.text.x  = element_blank(),axis.ticks.x  = element_blank(),
#                                        legend.position = 'left')
# p3=ggplot(tt_plot2[which(tt_plot2$response=='nonMHR'),],aes(y=cells,x=Abundance,fill=timepoint,label=pvals))+geom_bar(stat = 'identity',position = position_dodge(width = 0.8),color='black')+facet_grid(.~response)+
#   scale_fill_manual(values=c( "#AA4499","#DDCC77"),'Timepoint')+theme_classic(base_size = 15)+
#   theme(axis.text.x = element_text(angle = 90,hjust = 0.95,vjust = 0.5,color='black'),axis.text.y = element_text(color='black'))+xlab('% in all T cells')+ylab('')+
#   geom_text(size=8)
# p4=ggplot(tt_plot2[15:28,],aes(y=cells,fill=odds,x=-0.5))+geom_tile(color='grey30')+scale_fill_gradientn(colours = c("navy", "white", "firebrick"),limits = c(min(tt_plot2$odds[1:28]), max(tt_plot2$odds[1:28])),
#                                                                                                          values = rescale(c(min(tt_plot2$odds[1:28]),0, max(tt_plot2$odds[1:28]))),name='Log2\nodds\nratio')+
#   theme_classic(base_size = 15 )+theme(axis.title.x = element_blank(),axis.text.x  = element_blank(),axis.ticks.x  = element_blank(),legend.position = 'left')
# plot_grid(p2+ylab(''),
#           p1+theme(axis.title.y = element_blank(),axis.text.y  = element_blank(),axis.ticks.y  = element_blank(),axis.line.y = element_blank(),legend.position = 'none')+ylab(''),
#           p4+ylab(''),
#           p3+theme(axis.title.y = element_blank(),axis.text.y  = element_blank(),axis.ticks.y  = element_blank(),axis.line.y = element_blank(),legend.position = 'none')+ylab(''),
#           align = "h", ncol = 4, axis = "tb", rel_widths  = c(5, 10,5,10))

