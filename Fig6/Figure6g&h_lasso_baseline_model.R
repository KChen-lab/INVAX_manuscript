
########################################################Figure 6g&h########################

library(glmnet)
library(ggplot2)
library(readxl)


feature_data=read_xlsx('./Supplementary Tables 6-13_v3.xlsx',sheet = 'Supplementary Table 13',range = 'A2:FZ37')

#one-hot encoding for tumor stage variable
todummy=feature_data[,c(181:182)]
dummyresult <- fastDummies::dummy_cols(todummy, remove_first_dummy = F)
df=cbind(feature_data,dummyresult)
df$tobacco[which(df$tobacco==2)]=1
df <- df[!is.na(df$meta.Viable.change.percent), ]
#remove 33 because no baseline data, remove 26 because very few t cells at baseline
df <- df[which(df$INVAX!='INVAX026'&df$INVAX!='INVAX033'),]
df$Gender <- ifelse(df$Gender == 1, 0, 1)

y <- df$meta.Viable.change.percent

#get only baseline features
baseline.idx=grep("Baseline", colnames(df),value=T)
#include clinical features
baseline.idx=c(baseline.idx,c("Gender", "Age",'tobacco',grep('stage_',colnames(df),value = T)[2:8]))#do not include redundant stage variable

X.baseline <- df[, baseline.idx]
df.new <- cbind(y, X.baseline)
df.new <- df.new[complete.cases(df.new), ]
X.baseline <- as.matrix(df.new[, -1])
y <- df.new[, 1]

X.baseline <- scale(X.baseline)
X.baseline[, 2] <- df.new[, 3]

set.seed(1234)
fit <- cv.glmnet(x = X.baseline, y = y,
                 nlambda = 1e4,
                 type.measure = "mae",
                 standardize = FALSE)
fit$glmnet.fit$dev.ratio[which(fit$lambda == fit$lambda.min)]# 0.8921818


###Figure 6g scatter plot
y_pred <- predict(fit, newx = X.baseline)

df_plot=data.frame(Fitted=y_pred,real=y)
ggplot(df_plot,aes(x=lambda.1se,y=real))+geom_point()+theme_classic(base_size = 15)+ylab('Viable change %')+xlab('Predicted viable change %')+
  geom_abline(intercept = 0, slope = 1,color='#449979')+ylim(c(-100,-10))+xlim(c(-100,-10))+ggtitle(paste0('R-sqaured: ',format(fit$glmnet.fit$dev.ratio[which(fit$lambda == fit$lambda.min)] ,digits=4)))
ggsave('./Figure6h.scatter_lasso.pdf',width = 4,height = 3)


###Figure 6h bar plot
beta <- coef(fit, s = "lambda.min")
bb.nonzero <- beta@x
names(bb.nonzero) <- beta@Dimnames[[1]][beta@i + 1]
beta.inf <- coef(fit, x = X.baseline, y = y, s = "lambda.min", exact = TRUE)[-1]
names(beta.inf) <- beta@Dimnames[[1]][-1]
n=names(beta.inf)[beta.inf != 0]
n=str_replace(n,'.GSD','')
n[8]='CD16+ NK_Baseline_proportion'
results <- data.frame(Name = n,Coef=beta.inf[beta.inf != 0])
results$Name=factor(results$Name,level=results$Name[order(results$Coef,decreasing = TRUE)])
ggplot(results,aes(x=Coef,y=Name))+geom_bar(stat='identity',fill='#449979')+theme_classic(base_size = 12)+
  xlab('Lasso coefficients')+ylab('')+
  ggtitle(paste0('R-sqaured: ',format(fit$glmnet.fit$dev.ratio[which(fit$lambda == fit$lambda.min)] ,digits=4)))+theme(axis.text = element_text(colour = 'black'))
fit$glmnet.fit$dev.ratio[which(fit$lambda == fit$lambda.min)]+theme(axis.text = element_text(colour = 'black'))
ggsave('./barplot_coefficient.pdf',width = 10,height = 5)

