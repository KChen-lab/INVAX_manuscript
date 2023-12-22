######################################This code use GSDensity method to calculate gene set scores for CD8 T cells############
######################################We included gene set from Hallmark, KEGG, Reactome, and 3 defined gene expression signature###
######################################It need 

library(gsdensity)
library(supraHex)
library(RANN)
library(dnet)
library(anticlust)
library(multimode)
library(philentropy)
library(CelliD)
library(Seurat)
library(gridExtra)
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)
library(reshape2)
library(gsdensity)
library(supraHex)
library(RANN)
library(dnet)
library(anticlust)
library(multimode)
library(philentropy)
library(CelliD)
library(Seurat)
library(gridExtra)
library(msigdbr)
source("./GSDensity_code.R")


checkpoint=c('PDCD1','TIGIT','HAVCR2','CTLA4','LAG3','TOX')
cytotoxic=c('PRF1','GZMA','GZMB','GZMH','GNLY')
tissue_memory=c('ITGAE','ZNF683','ITGA1','CD69','CXCR6','CXCL13')

genesetlist=list(checkpoint,cytotoxic,tissue_memory)
genesetlist_names=c('checkpoint','cytotoxic','tissue_memor')
names(genesetlist) = genesetlist_names

msigdbr_df <- msigdbr(species = "human", category = "H")
pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

msigdbr_df <- msigdbr(species = 'Homo sapiens',category='C2',subcategory = 'CP:KEGG')
pathwaysK = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

msigdbr_df <- msigdbr(species = 'Homo sapiens',category='C2',subcategory = 'CP:REACTOME')
pathwaysR = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)


#if want to calculate pathway scores from public databases
big_list = c(genesetlist,pathwaysH,pathwaysK,pathwaysR)

#if want to calculate the checkpoint,cytotoxic and tissue_memory
big_list = c(genesetlist)


#read CD8 seurat object
data = readRDS("./CD8.rds")
data = subset(data,features = rownames(data)[-grep("^MT-",rownames(data))])


# compute cell/gene embeddings
ce <- compute.mca(object = data)

cells <- colnames(data)
el <- compute.nn.edges(coembed = ce, nn.use = 1000)

score_matrix = matrix(0,ncol(data),length(big_list))
for (i in 1:length(big_list)){
  cv <- run.rwr(el = el, gene_set = big_list[[i]], cells = cells)
  score_matrix[,i] = cv
}
#dev.off()
colnames(score_matrix) = names(big_list)
rownames(score_matrix) = cells
write.csv(score_matrix,paste0(path_result,"./GSDensity_score_CD8.csv"))

          
          