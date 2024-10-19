########################################################Extended Figure 10 ########################

library(glmnet)

df=read_xlsx('./Supplementary Tables 6-13_v3.xlsx',sheet = 'Supplementary Table 13',range = 'A2:FZ37')

df <- df[!is.na(df$meta.Viable.change.percent), ]
#remove 33 because no baseline data, remove 26 because very few t cells at baseline
df <- df[which(df$INVAX!='INVAX026'&df$INVAX!='INVAX033'),]#optional


y <- df$meta.Viable.change.percent
baseline.idx <- grepl("Baseline", colnames(df))
treatment.idx <- grepl("OnTreatment", colnames(df))
X.baseline <- df[, baseline.idx]
X.treatment <- df[, treatment.idx]

# Remove Baseline and OnTreatment strings
colnames(X.baseline) <- gsub("Baseline_", "", colnames(X.baseline))
colnames(X.baseline) <- gsub(".Baseline", "", colnames(X.baseline))
colnames(X.treatment) <- gsub("OnTreatment_", "", colnames(X.treatment))
colnames(X.treatment) <- gsub(".OnTreatment", "", colnames(X.treatment))

# Order columns of X.treatment according to X.baseline
X.treatment <- X.treatment[, match(colnames(X.treatment), colnames(X.baseline))]

# Create X.delta: X.treatment - X.baseline
X.delta <- X.treatment - X.baseline

No_of_expanded_clones <- df$Fraction_of_expanded_clones
No_of_contracted_clones <- df$Fraction_of_contracted_clones
df.new <- cbind(y,  X.delta, No_of_expanded_clones, No_of_contracted_clones)

df.new <- df.new[complete.cases(df.new), ]
colnames(df.new)

X.delta <- as.matrix(df.new[, -1])
y <- df.new[, 1]

# Scale X matrix
X.delta <- scale(X.delta)

set.seed(1234)
fit <- cv.glmnet(x = X.delta, y = y,
                 nlambda = 1e4,
                 type.measure = "mae",
                 standardize = FALSE)
fit$glmnet.fit$dev.ratio[which(fit$lambda == fit$lambda.min)] # 0.9823734
beta <- coef(fit, s = "lambda.min")
beta
beta.inf <- coef(fit, x = X.delta, y = y, s = "lambda.min", exact = TRUE)[-1]
lambda.inf <- fit$lambda.min * length(y)
names(beta.inf) <- beta@Dimnames[[1]][-1]

results <- data.frame(Name = names(beta.inf)[beta.inf != 0],
                      Coef = beta.inf[beta.inf != 0])


# Scatterplot
yhat <- predict(fit, s = "lambda.min", newx = X.delta)
df_plot=data.frame(Fitted=yhat,real=y)
ggplot(df_plot,aes(x=yhat,y=real))+geom_point()+theme_classic(base_size = 15)+ylab('Viable change %')+xlab('Predicted viable change %')+
  geom_abline(intercept = 0, slope = 1,color='#449979')+ylim(c(-100,-10))+xlim(c(-100,-10))+ggtitle(paste0('R-sqaured: ',format(fit$glmnet.fit$dev.ratio[which(fit$lambda == fit$lambda.min)] ,digits=4)))
ggsave('./scatterplot_delta_modeL.pdf',width = 4,height = 3)

# barplot for coefficients
results$Name <- factor(results$Name,level=results$Name[order(results$Coef,decreasing = TRUE)])
ggplot(results,aes(x = Coef, y = Name)) + 
  geom_bar(stat = 'identity', fill = '#449979') + 
  theme_classic(base_size = 15) +
  xlab('Lasso coefficients') + 
  ylab('') + 
  ggtitle(paste0('R-sqaured: ', format(fit$glmnet.fit$dev.ratio[which(fit$lambda == fit$lambda.min)], digits = 4)))
ggsave('./coeffecient_barplot_delta_model.pdf',width = 10,height = 5)
