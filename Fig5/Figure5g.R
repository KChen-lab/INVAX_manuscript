library(epitools)
library(ggplot2)
testProp <- function(metadata,timepoint,type){
  
  count_cluster.1=rev(table(metadata$CellSubset[which(metadata$Timepoint==timepoint&!is.na(metadata$vdj_barcode)&metadata$Response=='MHR')]))
  count_cluster.2=rev(table(metadata$CellSubset[which(metadata$Timepoint==timepoint&!is.na(metadata$vdj_barcode)&metadata$Response=='nonMHR')]))
  
  dat=metadata[which(metadata$Timepoint==timepoint&metadata$Expansion_annot==type&!is.na(metadata$Response)),]
  cellType=rev(levels(dat$CellSubset))
  
  stat.test <- data.frame(cellType = cellType)
  stat.test$group1 <- "MHR"
  stat.test$group2 <- "nonMHR"
  dat_ER=dat[which(dat$Response=='MHR'),]
  dat_R=dat[which(dat$Response=='nonMHR'),]
  stat.test$group1.size <- nrow(dat_ER)
  stat.test$group2.size <- nrow(dat_R)
  stat.test$group1_2.size <- nrow(dat_ER)+nrow(dat_R)
  stat.test$n1 <- NA 
  stat.test$n2 <- NA
  stat.test$prop1 <- NA
  stat.test$prop2 <- NA
  stat.test$prop1_all <- NA
  stat.test$prop2_all <- NA
  stat.test$prop1_cluster <- NA
  stat.test$prop2_cluster <- NA
  stat.test$prop1_cluster_all <- NA
  stat.test$prop2_cluster_all <- NA
  
  stat.test$statistic <- stat.test$df <- stat.test$p <- stat.test$odds <- stat.test$CI_d <- stat.test$CI_u <- NA
  stat.test$statistic_cluster <- stat.test$df_cluster <- stat.test$p_cluster <- stat.test$odds_cluster <- stat.test$cluster_CI_d <- stat.test$cluster_CI_u <- NA
  
  i <- 1
  for (cell.type in stat.test$cellType) {
    stat.test$n1[i] <- sum(dat_ER$CellSubset == cell.type)
    stat.test$n2[i] <- sum(dat_R$CellSubset == cell.type)
    stat.test$prop1[i] <- stat.test$n1[i] / nrow(dat_ER) *100
    stat.test$prop2[i] <- stat.test$n2[i] / nrow(dat_R) *100
    stat.test$prop1_all[i] <- stat.test$n1[i] / (nrow(dat_ER)+nrow(dat_R)) *100
    stat.test$prop2_all[i] <- stat.test$n2[i] / (nrow(dat_ER)+nrow(dat_R)) *100
    stat.test$prop1_cluster[i]=stat.test$n1[i] /count_cluster.1[cell.type] *100
    stat.test$prop2_cluster[i]=stat.test$n2[i] /count_cluster.2[cell.type] *100
    stat.test$prop1_cluster_all[i]=stat.test$n1[i] /(count_cluster.1[cell.type]+count_cluster.2[cell.type]) *100
    stat.test$prop2_cluster_all[i]=stat.test$n2[i] /(count_cluster.1[cell.type]+count_cluster.2[cell.type]) *100
    stat.test$n_patients1[i]=length(unique(dat_ER$patient[which(dat_ER$CellSubset == cell.type)]))
    stat.test$n_patients2[i]=length(unique(dat_R$patient[which(dat_R$CellSubset == cell.type)]))
    
    #test, among expanded/contracted clones, if more ER is found in cell.type                                                                                  
    # test <- prop.test(x = with(stat.test, c(n1[i], n2[i])), n = c(nrow(dat_ER), nrow(dat_R)))
    # stat.test$statistic[i] <- test$statistic
    # stat.test$df[i] <- test$parameter
    # stat.test$p[i] <- test$p.value
    test=oddsratio(matrix(c(stat.test$n1[i],nrow(dat_ER)-stat.test$n1[i],stat.test$n2[i],nrow(dat_R)-stat.test$n2[i]),nrow = 2),method = 'fisher')
    #stat.test$odds[i]<-log2(test$measure[2,1]/test$measure[1,1])
    test=fisher.test(matrix(c(stat.test$n1[i],nrow(dat_ER)-stat.test$n1[i],stat.test$n2[i],nrow(dat_R)-stat.test$n2[i]),nrow = 2))
    stat.test$p[i]=test$p.value
    stat.test$odds[i]=log2(test$estimate)
    stat.test$CI_d[i]=log2(test$conf.int[1])
    stat.test$CI_u[i]=log2(test$conf.int[2])
    
    i <- i + 1
  }
  
  stat.test$p.adj <- p.adjust(stat.test$p, method = "fdr")
  stat.test$p.adj_cluster <- p.adjust(stat.test$p_cluster, method = "fdr")
  
  for (i in 1:nrow(stat.test)){
    if (is.na(stat.test$p.adj[i])) {
      stat.test$p.adj.signif[i] <- "NA"
      # } else if (stat.test$p.adj[i] <= 1e-4) {
      #   stat.test$p.adj.signif[i] <- "****"
      # } else if (stat.test$p.adj[i] <= 0.001) {
      #   stat.test$p.adj.signif[i] <- "***"
      # } else if (stat.test$p.adj[i] <= 0.01) {
      #   stat.test$p.adj.signif[i] <- "**"
    } else if (stat.test$p.adj[i] <= 0.05) {
      stat.test$p.adj.signif[i] <- "*"
    } else {
      stat.test$p.adj.signif[i] <- "ns"
    }
  }
  stat.test <- tibble::as_tibble(stat.test)
  stat.test$cellType <- factor(stat.test$cellType, levels = rev(cellType))
  # stat.test$ymax <- pmax(stat.test$prop1, stat.test$prop2) + 0.005
  # stat.test$annotations <- sapply(1:nrow(stat.test), function(i) 
  #   paste(ifelse(stat.test$p.adj[i] <= 1e-3, format(stat.test$p.adj[i], scientific = TRUE, digits = 3), round(stat.test$p.adj[i], digits = 3)), "\n", stat.test$p.adj.signif[i]))
  return(stat.test)
}





#This code needs:
#scTCR-seq dataset includes cell barcode pseudoclonotype,clonotype (defined by CDR3A+CDR3B) and 
#whether clonotypes has expanded or not; the expansion test is done in Clonotyope_Expansion.R 
#cd8 meta data includes patient id, response,timepoint

data=readRDS('./CD8.rds')
wholevdj.dat=readRDS('./VDJ_data.rds')
#{r clonotype size umap }
meta=data@meta.data

#map scTCR data to cd8 data
ind1=match(intersect(wholevdj.dat$barcode,rownames(meta)),wholevdj.dat$barcode)
ind2=match(intersect(wholevdj.dat$barcode,rownames(meta)),rownames(meta))
identical(wholevdj.dat$barcode[ind1],rownames(meta)[ind2])
meta$Expansion_annot=NA
meta[ind2,'Expansion_annot']=wholevdj.dat[ind1,'Expansion_annot']
meta[ind2,'pseudoclonotype']=wholevdj.dat$pseudo_clonotype[ind1]
meta[ind2,'vdj_barcode']=wholevdj.dat$barcode[ind1]


type='Expanded_Sig'
cellType_marginal_test.b <- testProp(meta,1,type) # test difference in each cell type

cellType_marginal_test.onT <- testProp(meta,2,type) # test difference in each cell type


colors=c("#A6CEE3", "#4F96C4","royalblue" ,"#CC3311", "#EE3377", "#FB9A99", "#A7D78D","#69BB54","#5D9E43" ,"#117733" , "#DDCC77", "#D9A398", "#CAB2D6" ,"#9A77B8", "#6A3D9A"  )  

colnames(cellType_marginal_test.b)
cellType_marginal_test.b$cellType=factor(cellType_marginal_test.b$cellType,levels=levels(cellType_marginal_test.b$cellType))
df_plot=data.frame(celltype=c(as.character(cellType_marginal_test.b$cellType),as.character(cellType_marginal_test.onT$cellType)),odds=c(cellType_marginal_test.b$odds,cellType_marginal_test.onT$odds),
                   timepoint=c(rep('B',nrow(cellType_marginal_test.b)),rep('OnT',nrow(cellType_marginal_test.onT))),
                   pval=c(cellType_marginal_test.b$p,cellType_marginal_test.onT$p),
                   CI_b=c(cellType_marginal_test.b$CI_d,cellType_marginal_test.onT$CI_d),
                   CI_u=c(cellType_marginal_test.b$CI_u,cellType_marginal_test.onT$CI_u),
                   pval_adjust=c(cellType_marginal_test.b$p.adj,cellType_marginal_test.onT$p.adj))

df_plot$celltype=factor(df_plot$celltype,levels = rev(levels(cellType_marginal_test.b$cellType)))
ggplot(df_plot,aes(y=celltype,x=odds,color=celltype))+geom_boxplot()+
  geom_errorbar(aes(xmax=CI_u,xmin=CI_b),width=0.2)+geom_vline(xintercept = 0,linetype=2)+
  theme_minimal()+scale_color_manual(values=rev(colors))+facet_wrap(.~timepoint)+
  theme(legend.position = 'none',panel.background=element_rect(colour="black"))+ylab('')+xlab('Log2 odds ratio:MHR to nMHR')

############This code can be also used for extended figure 8a,d,h by changing the type name
