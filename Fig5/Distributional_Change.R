

wholevdj.dat=readRDS('~/Documents/ChenLab/Projects/HeadNeck/final.figures/Manuscript/Data and code/SeuratObjects for reviewers/VDJ_data.rds')
meta=readRDS('~/Documents/ChenLab/Projects/HeadNeck/final.figures/Manuscript/Data and code/SeuratObjects for reviewers/CD8_metadata.rds')


#map scTCR data to cd8 data
ind1=match(intersect(wholevdj.dat$barcode,rownames(meta)),wholevdj.dat$barcode)
ind2=match(intersect(wholevdj.dat$barcode,rownames(meta)),rownames(meta))
identical(wholevdj.dat$barcode[ind1],rownames(meta)[ind2])
meta$Expansion_annot=NA
meta[ind2,'Expansion_annot']=wholevdj.dat[ind1,'Expansion_annot']
meta[ind2,'pseudoclonotype']=wholevdj.dat$pseudo_clonotype[ind1]
meta[ind2,'vdj_barcode']=wholevdj.dat$barcode[ind1]


load("~/Documents/ChenLab/Projects/HeadNeck/Trajectory/heatmap_test_contract_sig.RData")
meta_large <- readRDS("~/Documents/ChenLab/Projects/HeadNeck/final.figures/Manuscript/Data and code/cd8_meta_newexpansion.rds")
meta <- meta[, c("pseudoclonotype", "CellSubset", "Timepoint", "Expansion_annot")]
meta <- meta[complete.cases(meta), ]
meta <- subset(meta, Expansion_annot == "Expanded_Sig")
colnames(meta) <- c("clonotype", "cell.type", "timepoint", "expansion")
meta$clonotype <- as.factor(meta$clonotype)
meta$timepoint <- as.factor(meta$timepoint)

data_pre <- subset(meta, timepoint == 1)
data_post <- subset(meta, timepoint == 2)

d_pre <- table(data_pre$clonotype, data_pre$cell.type)
d_post <- table(data_post$clonotype, data_post$cell.type)
all.equal(rownames(d_pre), rownames(d_post))
all.equal(colnames(d_pre), colnames(d_post))

pval_cell_dist <- rep(NA, nrow(d_post))
for (i in 1:nrow(d_post)) {
  x <- rbind(d_pre[i, ], d_post[i, ])
  test <- fisher.test(x, workspace = 10e7)
  pval_cell_dist[i] <- test$p.value
}
pval.adj_cell_dist <- p.adjust(pval_cell_dist, method = "fdr")
sum(pval.adj_cell_dist <= 0.05, na.rm = TRUE)
names(pval_cell_dist) <- names(pval.adj_cell_dist) <- rownames(d_post)


pval_celltype_dist <- data.frame(pval = pval_cell_dist, pval.adj = pval.adj_cell_dist)

