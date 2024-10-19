

#############################################Figure 3i&h patients with progression or not ####################

library(ggplot2)
library(readxl)
library(reshape2)
library(ggpubr)
library(stringr)
library(rstatix)
feature_data=read_xlsx('./Supplementary Tables 6-13_v3.xlsx',sheet='Supplementary Table 14',
                       range='A2:FZ37')


#get only CD4, CD8 subset and CD16 NK data.
f=colnames(feature_data)[c(44:101,108,179)]

df=feature_data[,c(colnames(feature_data)[c(1,2,3)],f)]
colnames(df)

long_df <- melt(df, id.vars = colnames(feature_data[c(1,2,3)]), 
                variable.name = "cellsubset", value.name = "proportions")


long_df=long_df[which(!is.na(long_df$Progress)),]
#baseline data for patient 26 was removed due to small number of lymphocytes.
todel=which(long_df$INVAX=='26'&str_detect(long_df$cellsubset,'Baseline'))
long_df=long_df[-todel,]
long_df$proportions=long_df$proportions*100
long_df=long_df[which(!is.na(long_df$proportions)),]


ggplot(long_df[which(long_df$cellsubset=='CD16NK_proportion.Baseline'|long_df$cellsubset=='Treg.Effector_OnTreatment_proportion'),],
       aes(x=Progress,y=proportions,fill=Progress))+
  geom_boxplot()+stat_compare_means(label='p.format',size=5)+facet_wrap(.~cellsubset,scales = 'free')+geom_boxplot()+
  geom_point()+theme_classic(base_size = 10)+xlab('')+scale_fill_manual(values=c('forestgreen','orange'))+ylab('% Abundance')+
  theme(legend.position = 'none',strip.background = element_blank(), axis.text.x = element_text(size=20),strip.text = element_text(size=12),
        axis.text.y = element_text(size=20) ,axis.title.y = element_text(size=21))


long_df.sub=long_df[which(long_df$cellsubset=='CD16NK_proportion.Baseline'|long_df$cellsubset=='Treg.Effector_OnTreatment_proportion'),]
stat.test=wilcox_test(proportions~Progress,data=long_df.sub[which(long_df.sub$cellsubset=='CD16NK_proportion.Baseline'),])
stat.test=wilcox_test(proportions~Progress,data=long_df.sub[which(long_df.sub$cellsubset=='Treg.Effector_OnTreatment_proportion'),])




#############################################This code can also be used for supplementary figure 6 ####################

#plots for all CD4, CD8 subsets at baseline
long_df.baseline=long_df[which(str_detect(long_df$cellsubset,'Baseline')),]
long_df.baseline$cellsubset=as.character(long_df.baseline$cellsubset)
long_df.baseline$cellsubset=str_replace(long_df.baseline$cellsubset,'_Baseline_proportion','')
long_df.baseline$cellsubset=str_replace(long_df.baseline$cellsubset,'_proportion.Baseline','')

length(unique(long_df.baseline$INVAX))#28 patients;if only focus on patients with histological response data, 24 patients.
long_df.baseline$cellsubset=factor(long_df.baseline$cellsubset,levels = unique(long_df.baseline$cellsubset))
ggplot(long_df.baseline, aes(x=Progress,y=proportions,fill=Progress))+
  geom_boxplot()+stat_compare_means(label='p.format',size=5)+facet_wrap(.~cellsubset,ncol =10)+geom_boxplot()+
  geom_point()+theme_classic(base_size = 10)+xlab('')+scale_fill_manual(values=c('forestgreen','orange'))+ylab('% Abundance')+
  theme(legend.position = 'none',strip.background = element_blank(), axis.text.x = element_blank(),strip.text = element_text(size=12),
        axis.text.y = element_text(size=20) ,axis.title.y = element_text(size=21))

#plots for all CD4, CD8 subsets on treatment (on induction)
long_df.ont=long_df[which(str_detect(long_df$cellsubset,'OnTreatment')),]
length(unique(long_df.ont$INVAX))#30 patients; if only focus on patients with histological response data, 26 patients.
long_df.ont$cellsubset=as.character(long_df.ont$cellsubset)
long_df.ont$cellsubset=str_replace(long_df.ont$cellsubset,'_OnTreatment_proportion','')
long_df.ont$cellsubset=str_replace(long_df.ont$cellsubset,'_proportion.OnTreatment','')

long_df.ont$cellsubset=factor(long_df.ont$cellsubset,levels = unique(long_df.ont$cellsubset))
ggplot(long_df.ont, aes(x=Progress,y=proportions,fill=Progress))+
  geom_boxplot()+stat_compare_means(label='p.format',size=5)+facet_wrap(.~cellsubset,ncol =10)+geom_boxplot()+
  geom_point()+theme_classic(base_size = 10)+xlab('')+scale_fill_manual(values=c('forestgreen','orange'))+ylab('% Abundance')+
  theme(legend.position = 'none',strip.background = element_blank(), axis.text.x = element_blank(),strip.text = element_text(size=12),
        axis.text.y = element_text(size=20) ,axis.title.y = element_text(size=21))


  
#statistical test
perform_wilcox_test <- function(cellsubset_value, data) {
  filtered_data <- data[data$cellsubset == cellsubset_value, ]
  test_result <- wilcox_test(proportions ~ Progress, data = filtered_data)
  return(test_result)
}

cellsubset_levels <- unique(long_df.baseline$cellsubset)
test_results <- lapply(cellsubset_levels, perform_wilcox_test, data = long_df.baseline)
combined_results <- do.call(rbind, test_results)
combined_results=as.data.frame(combined_results)
combined_results$p.adj <- p.adjust(combined_results$p, method = "fdr")
rownames(combined_results)=as.character(cellsubset_levels)


cellsubset_levels <- unique(long_df.ont$cellsubset)
test_results <- lapply(cellsubset_levels, perform_wilcox_test, data = long_df.ont)
combined_results <- do.call(rbind, test_results)
combined_results=as.data.frame(combined_results)
combined_results$p.adj <- p.adjust(combined_results$p, method = "fdr")
rownames(combined_results)=as.character(cellsubset_levels)
