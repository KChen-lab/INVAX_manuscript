########################################Figure 3c CTLA4-expression and depletion with statistical test############

library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(readxl)
library(ggplot2)
source("./tradeSeqFunc.R")
set.seed(123)


data=readRDS('./CD4.rds') #This code needs CD4 seurat object

ctla4=FetchData(subset(data,idents = 'Treg.Effector',subset=PATH.Response.new=='ER'|PATH.Response.new=='R'),vars = 'CTLA4')
ctla4.ER.2=FetchData(subset(data,idents = 'Treg.Effector',subset=PATH.Response.new=='ER'&Timepoint==2),vars = 'CTLA4')
ctla4.ER.1=FetchData(subset(data,idents = 'Treg.Effector',subset=PATH.Response.new=='ER'&Timepoint==1),vars = 'CTLA4')
ctla4.R.2=FetchData(subset(data,idents = 'Treg.Effector',subset=PATH.Response.new=='R'&Timepoint==2),vars = 'CTLA4')
ctla4.R.1=FetchData(subset(data,idents = 'Treg.Effector',subset=PATH.Response.new=='R'&Timepoint==1),vars = 'CTLA4')


q=seq(0,100,1)

r1=c()
r2=c()
r3=c()
r4=c()

for(i in seq(0,99,1)){
th=quantile(ctla4$CTLA4,probs=i/100)
r1=c(r1,length(which(ctla4.ER.2$CTLA4>th))/nrow(ctla4.ER.2))
r2=c(r2,length(which(ctla4.ER.1$CTLA4>th))/nrow(ctla4.ER.1))
r3=c(r3,length(which(ctla4.R.2$CTLA4>th))/nrow(ctla4.R.2))
r4=c(r4,length(which(ctla4.R.1$CTLA4>th))/nrow(ctla4.R.1))
}

dat=data.frame(qua=seq(0,99,1),ER=r1-r2,R=r3-r4)
dat=data.frame(val=c((r1-r2)/r2*100,(r3-r4)/r4*100),qua=c(seq(0,99,1),seq(0,99,1)),
               response=c(rep('MHR',100),rep('nMHR',100)))


countMatrix <- matrix(dat$val, nrow = 1)
rownames(countMatrix) <- "n1"
colnames(countMatrix) <- paste0("c", 1:200)
conditions <- as.factor(dat$response)
pseudotime <- as.matrix(dat$qua)
cellWeights <- matrix(1, nrow = nrow(dat))
rownames(pseudotime) <- rownames(cellWeights) <- colnames(countMatrix)


# Wald test
sce <- fitGAM(as.matrix(countMatrix), pseudotime = pseudotime,
              nknots = 20,
              cellWeights = cellWeights,
              conditions = conditions, 
              family = "gaussian",
              offset = rep(0, ncol(dat)) ,
              sce = TRUE
)
(res <- conditionTest(sce$sc))



# Plot fitted models & CIs
model <- sce$gamList$n1
 p_obj <- plot(model, shade = TRUE, seWithMean = TRUE, residuals = TRUE, 
              pch = 16, cex = 0.8)
# p_obj
p1 <- as.data.frame(p_obj[[1]][c("x", "se", "fit", "raw", "p.resid")])[1:100, ]
p2 <- as.data.frame(p_obj[[2]][c("x", "se", "fit", "raw", "p.resid")])[101:200, ]

gam_df <- rbind(p1, p2)
rownames(gam_df) <- NULL
gam_df$response <- rep(c("MHR", "nMHR"), each = 100)


ggplot(gam_df, aes(x = x, group = response, color = response)) +geom_hline(yintercept = 0)+geom_line(aes(x = x, y = fit)) +
  geom_point(aes(x = raw, y = p.resid)) +
  geom_ribbon(aes(ymin = fit - 2*se, ymax = fit + 2*se, y = NULL),
              alpha = 0.3) +
  scale_color_manual(values =c('royalblue','lightblue'))+
  labs(x = "Percentile of CTLA-4 expression level of all effector Treg cells",
       y = "% Reduciton of CTLA-4 expressing cell percentages")+theme_bw()

