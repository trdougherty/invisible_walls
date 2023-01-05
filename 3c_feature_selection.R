library(glmnet)
library(caret)
library(car)
library(readr)
library(dplyr)
library(fastDummies)
library(corrplot)
# library(Hmisc)
library(mclust)
library(ClusterR)
library("factoextra")
library(ggcorrplot)
library(bigmemory)
library(tibble)
library(tidyr)
library(scales)
library(xtable)
library(texreg)
library(stargazer)
library(Rfast)
library(gplots)

library("RColorBrewer")
library("FactoMineR")


library(pacman)

pacman::p_install_gh("rstudio/d3heatmap")
pacman::p_load(webshot)

library(ggfortify)


# library(devtools)
# install_github("vqv/ggbiplot")
# library(ggbiplot)

library(doParallel)
registerDoParallel(cores=8)

source("3c_helper_functions.R")
set.seed(42)

# # exploring stats for the paper a bit
# footprints <- read_csv("../building_p.csv")
# quantile(footprints$HEIGHTROOF * 0.3048, 0.99, na.rm=TRUE)

training <- read_csv("3a/training.csv")
validating <- read_csv("3a/validating.csv")
testing <- read_csv("3a/testing.csv")

full_data <- rbind(training, validating, testing)
full_data

h <- hist(full_data$electricity_mwh, breaks=10)
text(h$mids,h$counts,labels=h$counts, adj=c(0.5, -0.5))

h <- hist(full_data$gas_mwh, breaks=10)
text(h$mids,h$counts,labels=h$counts, adj=c(0.5, -0.5))

mean(full_data$gas_mwh)

# Looking into scaled and normalized data to study variable effect size
###############
###############
###############

################# Now into the preprocessing for the regression
prepped_data <- prep_data(full_data, scale=FALSE)

electric_norm <- prepped_data[1]
gas_norm <- prepped_data[2]
training_prepped <- prepped_data[3]

### ELECTRIC DATA
train_electric <- data.frame(training_prepped)
train_electric$electric_norm <- as.numeric(unlist(electric_norm))

# want to drop the values but also want to keep track of the unique IDs we're losing
keep_idx_electric <- (train_electric$electric_norm == -Inf) | (train_electric$electric_norm == Inf) | !complete.cases(train_electric)
train_electric <- train_electric[!keep_idx_electric,]

electric_idx <- full_data[!keep_idx_electric,c("bbl","bin","year","month")]
electric_idx

train_electric
training_prepped

small_electric <- train_electric[1:10000,]

train_ex <- train_electric[, -which(names(train_electric) =="electric_norm")]
train_ex <- as.matrix(sapply(train_ex, as.numeric))

env_clustering_terms_pca <- PCA(train_ex, graph=FALSE)
var <- get_pca_var(env_clustering_terms_pca)

tiff(
  "/Users/thomas/Desktop/Work/Research/uil/postquals/thermal/data/general_pca.tiff",
  units="cm", 
  width=21, 
  height=21, 
  res=300, 
  pointsize=6
)

fviz_pca_var(
  env_clustering_terms_pca,
  axes = c(1, 2),
  repel = TRUE,
  alpha.var = "contrib"
)
dev.off()

dup_cor <- cor(train_ex)
ggcorrplot::ggcorrplot(dup_cor)
# 
col<- colorRampPalette(c("white", "black"))(10)


cormat <- round(dup_cor,2)

# Get upper triangle of the correlation matrix
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)]<- NA
  return(cormat)
}

lower_tri <- get_lower_tri(cormat)

library(reshape2)
melted_cormat <- melt(lower_tri, na.rm = TRUE)

# Heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()


heatmap(
  x = dup_cor,
  col = col,
  trace="none",
  symm = TRUE,
)

tiff("/Users/thomas/Desktop/Work/Research/uil/postquals/thermal/data/corrplot.tiff", units="cm", width=18, height=18, res=600, pointsize=10)
heatmap.2(
  dup_cor,
  col=brewer.pal(9,"YlGn"),
  dendrogram='row', 
  Rowv=TRUE, 
  Colv=TRUE, 
  trace="none"
)
dev.off()



heatmap.2(
  dup_cor, 
  symm = TRUE,
  col=bluered, 
  scale="column",
  tracecol="#303030"
)


d3heatmap(
  dup_cor, 
  dendrogram = 'none',
  key.title = "Legend",
  scale = "row",
  colors = "Spectral",
  type="upper"
)

# correlated_idx <- findCorrelation(dup_cor, cutoff = 0.9)
# dup_cor[,correlated_idx]
# 
# cor_matrix[,correlated_idx]
# 
# colnames(train_ex[,c(correlated_idx)])[1:15]
# 
# ggcorrplot::ggcorrplot(dup_cor[,colnames(train_ex[,c(correlated_idx)])[1:14]][70:90,])
# 

# want to remove the columns which have a high level of correlation
correlated_features <- c(
  # Term 1, represented by B1
  # "B11",
  "B8A",
  "B7",
  "B8",
  "B6",
  "B2",
  "B12",
  "B9",
  "B3",
  "B4",
  "B5",
  # Term 2, represented by HDD
  "SPFH",
  "DPT",
  "TMP",
  # Term 3, represented by WIND
  "GUST",
  "UGRD",
  "VGRD",
  "WDIR",
  "PRES",
  "VIS",
  # Term 4, represented by elevation
  "groundelev",
  "HGT",
  # Term 5, represented by assesstot
  "assessland",
  "exempttot",
  # Term 6, represented by TCDC
  # "TCDC",
  "cf_cvg"
)

# select(data.frame(train_ex), correlated_features)
uncorrelated <- train_electric[,!(colnames(train_electric) %in% correlated_features)]

# stargazer(uncorrelated, out = "latex_code/data_summary.txt")

linear_model_electric_mc <- lm(electric_norm~., uncorrelated)
linear_vif_uncorr <- car::vif(linear_model_electric_mc)
barplot(linear_vif_uncorr, main = "VIF Values", horiz = TRUE, col = "indianred")

# dup_cor_uc <- cor(uncorrelated)
# dup_cor_uc[,"assesstot"]
# heatmap(x = dup_cor_uc, col = col, symm = TRUE)

vif_uncorr <- data.frame(linear_vif_uncorr[order(linear_vif_uncorr)])
vif_uncorr

# print(xtable(round(vif_uncorr, digits=0), type = "latex"), file = "latex_code/vif_table_electric_uncorr.tex")

# now I want to try and actually interpret the regression
# summary(linear_model_electric_mc)

# reg1 <- lm(electric_norm~hdd,data=small_electric)
# with(small_electric,plot(hdd, electric_norm))
# abline(reg1, col="red")

# print(texreg(linear_model_electric_mc), file="latex_code/electric_regression_normalized.tex")

###
### GAS DATA, also normed by std.dev
train_gas <- data.frame(training_prepped)
train_gas$gas_norm <- as.numeric(unlist(gas_norm))
# want to drop the values but also want to keep track of the unique IDs we're losing

keep_idx_gas <- (train_gas$gas_norm == -Inf) | (train_gas$gas_norm == Inf) | !complete.cases(train_gas)
train_gas <- train_gas[!keep_idx_gas,]

train_ex_gas <- train_gas[, -which(names(train_gas) =="gas_norm")]
train_ex_gas <- as.matrix(sapply(train_ex_gas, as.numeric))

gas_idx <- full_data[!keep_idx_gas,]
gas_idx

# now to remove the aliased terms
asiased_features_gas <- c(
  "TMP"
)

# select(data.frame(train_ex), correlated_features)
train_gas <- train_gas[,!(colnames(train_gas) %in% asiased_features_gas)]
linear_model_gas_mc <- lm(gas_norm~., train_gas)

linear_vif_gas <- car::vif(linear_model_gas_mc)
barplot(linear_vif_gas, main = "VIF Values", horiz = TRUE, col = "black")

vif_gas_uncorr <- data.frame(linear_vif_gas[order(linear_vif_gas)])
vif_gas_uncorr

# select(data.frame(train_ex), correlated_features)
train_gas_uncorr <- train_gas[,!(colnames(train_gas) %in% correlated_features)]
# stargazer(train_gas_uncorr, out = "latex_code/data_summary_gas_normalized.txt")

linear_model_gas_mc <- lm(gas_norm~., train_gas_uncorr)
linear_vif_gas_uncorr <- car::vif(linear_model_gas_mc)
barplot(linear_vif_gas_uncorr, main = "VIF Values", horiz = TRUE, col = "black")


# print(xtable(round(vif_uncorr, digits=0), type = "latex"), file = "latex_code/vif_table_uncorr.tex")

# # now I want to try and actually interpret the regression
# summary(linear_model_gas_mc)
# 
# print(xtable(round(vif_uncorr, digits=0), type = "latex"), file = "latex_code/vif_table_uncorr_normalized_gas.tex")
# stargazer(linear_model_gas_mc, out ="latex_code/gas_regression_normalized.tex")

data.frame(format((exp(linear_model_electric_mc$coefficients) - 1) * 100, scientific = F, digits=2))

stargazer(
  linear_model_electric_mc, 
  linear_model_gas_mc,
  dep.var.labels=c("log( Electric / Area )","log( Gas / Area )"),
  title="Linear Model - Mean Centered",  
  ci=FALSE,
  single.row=TRUE,
  type="text"
)

capture.output(stargazer(
  linear_model_electric_mc, 
  linear_model_gas_mc,
  dep.var.labels=c("log( Electric / Area )","log( Gas / Area )"),
  title="Linear Model - Mean Centered",  
  ci=FALSE,
  single.row=TRUE,
  out= "latex_code/electric_gas_mean_centered.tex"
))


#################
#################
################# Switching to non-normalized std dev data
prepped_data <- prep_data(full_data, scale=TRUE)

electric_norm <- prepped_data[1]
gas_norm <- prepped_data[2]
training_prepped <- prepped_data[3]

### ELECTRIC DATA
train_electric <- data.frame(training_prepped)
train_electric$electric_norm <- as.numeric(unlist(electric_norm))

# want to drop the values but also want to keep track of the unique IDs we're losing
keep_idx_electric <- (train_electric$electric_norm == -Inf) | (train_electric$electric_norm == Inf) | !complete.cases(train_electric)
train_electric <- train_electric[!keep_idx_electric,]

train_ex <- train_electric[, -which(names(train_electric) =="electric_norm")]
train_ex <- as.matrix(sapply(train_ex, as.numeric))

# dup_cor <- cor(train_ex)
# heatmap(x = dup_cor, col = col, symm = TRUE)

# select(data.frame(train_ex), correlated_features)
uncorrelated <- train_electric[,!(colnames(train_electric) %in% correlated_features)]
length(colnames(uncorrelated))
capture.output(stargazer(uncorrelated, out = "latex_code/data_summary.txt"))

uncorrelated

linear_model_electric_norm <- lm(electric_norm~., uncorrelated)
linear_vif_uncorr <- car::vif(linear_model_electric_norm)
barplot(linear_vif_uncorr, main = "VIF Values", horiz = TRUE, col = "indianred")

# dup_cor_uc <- cor(uncorrelated)
# heatmap(x = dup_cor_uc, col = col, symm = TRUE)

vif_uncorr <- data.frame(linear_vif_uncorr[order(linear_vif_uncorr)])

# 
# # now I want to try and actually interpret the regression
# summary(linear_model_electric_norm)
# 
# print(xtable(round(vif_uncorr, digits=0), type = "latex"), file = "latex_code/vif_table_uncorr_electric.tex")
# print(texreg(linear_model_electric_norm), file="latex_code/electric_regression.tex")


### GAS DATA
train_gas <- data.frame(training_prepped)
train_gas$gas_norm <- as.numeric(unlist(gas_norm))

# want to drop the values but also want to keep track of the unique IDs we're losing
keep_idx_gas <- (train_gas$gas_norm == -Inf) | (train_gas$gas_norm == Inf) | !complete.cases(train_gas)
train_gas <- train_gas[!keep_idx_gas,]

train_ex_gas <- train_gas[, -which(names(train_gas) =="gas_norm")]
train_ex_gas_t <- as.matrix(sapply(train_ex_gas, as.numeric))
train_ex_gas_t

# now to remove the aliased terms
asiased_features_gas <- c(
  "TMP"
)

# select(data.frame(train_ex), correlated_features)
train_gas <- train_gas[,!(colnames(train_gas) %in% asiased_features_gas)]
linear_model_gas_norm <- lm(gas_norm~., train_gas)
linear_vif_gas <- car::vif(linear_model_gas_norm)
barplot(linear_vif_gas, main = "VIF Values", horiz = TRUE, col = "black")

vif_gas_uncorr <- data.frame(linear_vif_gas[order(linear_vif_gas)])

# select(data.frame(train_ex), correlated_features)
train_gas_uncorr <- train_gas[,!(colnames(train_gas) %in% correlated_features)]
stargazer(train_gas_uncorr, out = "latex_code/data_summary_gas_norm.txt")

linear_model_gas_norm <- lm(gas_norm~., train_gas_uncorr)
linear_vif_gas_uncorr <- car::vif(linear_model_gas_norm)
barplot(linear_vif_gas_uncorr, main = "VIF Values", horiz = TRUE, col = "black")


# print(xtable(round(vif_uncorr, digits=0), type = "latex"), file = "latex_code/vif_table_uncorr.tex")

# now I want to try and actually interpret the regression
# summary(linear_model_gas_norm)
# 
# print(xtable(round(vif_uncorr, digits=0), type = "latex"), file = "latex_code/vif_table_gas_norm.tex")
# print(texreg(linear_model_gas_norm), file="latex_code/gas_regression_norm.tex")

stargazer(
  linear_model_electric_norm, 
  linear_model_gas_norm, 
  title="Linear Model - Normalized",  
  ci=TRUE,
  dep.var.labels=c("log( Electric / Area )","log( Gas / Area )"),
  single.row=TRUE,
  model.numbers = FALSE,
  multicolumn = TRUE,
  type="text"
)

stargazer(
  linear_model_electric_norm, 
  linear_model_gas_norm, 
  title="Linear Model - Normalized",  
  ci=TRUE,
  dep.var.labels=c("log( Electric / Area )","log( Gas / Area )"),
  single.row=TRUE,
  model.numbers = FALSE,
  multicolumn = TRUE,
  out= "latex_code/electric_gas_normalized.tex"
)

###################
###################
###################
uncorrelated <- train_electric[,!(colnames(train_electric) %in% correlated_features)]
colnames(uncorrelated)

data_sources <- c(
  "VIIRS",
  "Sentinel2-1A",
  "Sentinel2-1A",
  "Sentinel2-1A",
  "NOAA",
  "NOAA",
  "NOAA",
  "NOAA",
  "NOAA",
  "NASA JPL",
  "PLUTO",
  "PLUTO"
)

uncorr_names <- c(colnames(uncorrelated)[1:length(colnames(uncorrelated))-1])

fortable <- lm(electric_norm~., uncorrelated)
fortable_vif <- car::vif(fortable)
barplot(fortable_vif, main = "VIF Values", horiz = TRUE, col = "black")


final_datatypes <- data.frame(uncorr_names, data_sources, fortable_vif)

final_datatypes
rownames(final_datatypes) <- 1:nrow(final_datatypes)
colnames(final_datatypes) <- c("Data Description", "Data Source", "VIF")
xtable(final_datatypes)




##### ELECTRIC NORM
train_ex_intermediate <- as.matrix(data.frame(train_ex)[,!(colnames(train_ex) %in% correlated_features)])
train_electric_df <- data.frame(train_ex_intermediate)
train_electric_df

train_electric_df$electric_norm <- train_electric$electric_norm

train_lm_electric <- lm(electric_norm~.,train_electric_df)
elastic_coef_electric <- coef(train_lm_electric)
# elastic_coef_electric
# 

# how we can predict with the linear model
elastic_coef_terms <- elastic_coef_electric[2:length(elastic_coef_electric)]

train_env_scaled_electric <- t(t(train_ex_intermediate) * elastic_coef_terms)
train_env_scaled_electric
# 
# 
predict(train_lm_electric, newdata = data.frame(train_ex_intermediate))

# train_env_vars <- train_ex[,2:train_ex_idxend]
# env_coef <- elastic_coef_electric[3:end_idx,]

env_coef <- colnames(train_ex)[1:31]
# 
# train_env_vars <- train_ex[,2:73]
# env_coef <- elastic_coef[3:74,1]
# 
# data.frame(elastic_coef_terms[order(elastic_coef_terms),drop=FALSE])
# 
# train_env_scaled <- t(t(train_env_vars) * env_coef)
# train_env_scaled
# 
train_env_scaled_electric
train_env_scaled_electric

# this comes without the intercept term
train_env_predictions <- data.frame(rowSums(train_env_scaled_electric))
train_env_predictions
colnames(train_env_predictions) <- c("predictions")
train_env_predictions

lasso_dropcols <- names(env_coef[env_coef == 0])
env_clustering_terms <- train_env_scaled_electric[,!(colnames(train_env_scaled_electric) %in% lasso_dropcols)]

exp(rowSums(env_clustering_terms)) - 1

env_clustering_terms_pca <- prcomp(env_clustering_terms, scale=FALSE)
var <- get_pca_var(env_clustering_terms_pca)

tiff("/Users/thomas/Desktop/Work/Research/uil/postquals/thermal/data/electric_pca.tiff", units="cm", width=15, height=15, res=300, pointsize=6)
fviz_pca_var(
  env_clustering_terms_pca,
  axes = c(1, 2),
  repel = TRUE,
  alpha.var = "contrib"
)
dev.off()


colnames(env_clustering_terms) <- paste(colnames(env_clustering_terms),"reg",sep="_")

# we want to now just strip off the environmental variables, which this should do
env_terms_only <- env_clustering_terms[,2:(ncol(env_clustering_terms)-3)]
colnames(env_terms_only)

env_terms_only

train_envonly_predictions <- data.frame(rowSums(env_terms_only))
colnames(train_envonly_predictions) <- c("env_predictions")
train_envonly_predictions
# 
group_3 <- kmeans(env_terms_only, centers=3)
group_5 <- kmeans(env_terms_only, centers=5)
group_10 <- kmeans(env_terms_only, centers=10)
# group_15 <- kmeans(env_terms_only, centers=15)
# group_20 <- kmeans(env_terms_only, centers=20)
# 
gmm_groups_electric <- Mclust(env_terms_only, G=5)
gmm_groups_electric$parameters$mean

library(d3heatmap)
d3heatmap(
  gmm_groups_electric$parameters$mean, 
  dendrogram = 'none',
  key.title = "Legend",
  print.values = T,
  scale = "column",
  colors = "RdYlBu"
)

heatmap(
  gmm_groups_electric$parameters$mean,
  dendrogram='none', 
  Rowv=FALSE, 
  Colv=FALSE,
  trace='none'
)

100 * (exp(colSums(gmm_groups_electric$parameters$mean)) - 1)

heatmap(gmm_groups_electric$parameters$mean)

tt <- data.frame(gmm_groups_electric$parameters[2])
tt

exp(sum(data.frame(gmm_groups_electric$parameters[2])[,"mean.1"]))
exp(sum(data.frame(gmm_groups_electric$parameters[2])[,"mean.2"]))


tt['diff'] <- abs(tt[,"mean.1"] - tt[,"mean.2"])

tt
# 
# summary(gmm_groups)
# plot(gmm_groups)
jpeg(
  "/Users/thomas/Desktop/Work/Research/uil/postquals/thermal/data/electric_gmm.jpeg",
  height=500,
  width=500,
  quality=300
)
fviz_cluster(
  gmm_groups_electric,
  axes = c(1, 2),
  data=env_terms_only,
  geom="point",
  ggtheme = theme_light(),
  pointsize=0.4,
  dpi=300
)
dev.off()

# train_returning <- cbind(
#   electric_idx,
#   gmm_groups$classification, 
#   env_terms_only, 
#   train_env_predictions, 
#   train_envonly_predictions
# )
# 
# write.csv(train_returning, "3c/training_environmental_gmm.csv", row.names = FALSE)

g5_cluster_electric <- as.matrix(gmm_groups_electric$classification)
colnames(g5_cluster_electric) <- c("cluster")


percent_deviation_electric <- 100* (exp(train_envonly_predictions) - 1)
colnames(percent_deviation_electric) <- c("percent_deviation_env")


d3heatmap(
  gmm_groups_electric$parameters$mean, 
  dendrogram = 'none',
  key.title = "Legend",
  print.values = T,
  scale = "row",
  colors = "Spectral"
)

electric_centroids <- data.frame(100 * (exp(colSums(gmm_groups_electric$parameters$mean)) - 1))
colnames(electric_centroids) <- c("Predicted Environmental Impact")
xtable(electric_centroids)

env_terms_only
train_returning_gmm <- cbind(
  electric_idx, 
  g5_cluster_electric, 
  env_terms_only,
  train_env_predictions, 
  train_envonly_predictions, 
  percent_deviation_electric
)

train_returning_gmm
write.csv(train_returning_gmm, "3c/training_environmental_electric_gmm.csv", row.names = FALSE)

train_returning_km <- cbind(
  electric_idx,
  group_5$cluster, 
  env_terms_only, 
  train_env_predictions, 
  train_envonly_predictions
)
write.csv(train_returning_km, "3c/training_environmental_km5.csv", row.names = FALSE)

##### GAS NORM
train_ex_gas_t
train_ex_intermediate_gas <- as.matrix(data.frame(train_ex_gas_t)[,!(colnames(train_ex_gas_t) %in% correlated_features)])
train_gas_df <- data.frame(train_ex_intermediate_gas)

train_gas_df$gas_norm <- train_gas$gas_norm

train_lm_gas <- lm(gas_norm~.,train_gas_df)
train_lm_gas

# train_elastic_gas = cv.glmnet(
#   train_ex_gas_t,
#   train_gas$gas_norm,
#   type.measure = "mse",
#   nfolds = 5,
#   parallel = TRUE
# )
# 
elastic_coef_gas <- coef(train_lm_gas)
elastic_coef_gas_terms <- elastic_coef_gas[2:length(elastic_coef_gas)]
elastic_coef_gas_terms
# 
# how we can predict with the linear model
predict(train_lm_gas, newdata = data.frame(train_ex_gas_t))

train_ex_intermediate_gas

train_env_scaled_gas <- t(t(train_ex_intermediate_gas) * elastic_coef_gas_terms)
train_env_scaled_gas

# train_env_vars <- train_ex[,2:train_ex_idxend]
# env_coef <- elastic_coef[3:end_idx,]

# colnames(train_ex)[2:73]
# 
# train_env_vars <- train_ex[,2:73]
# env_coef <- elastic_coef[3:74,1]
# 
# data.frame(elastic_coef_terms[order(elastic_coef_terms),drop=FALSE])
# 
# train_env_scaled <- t(t(train_env_vars) * env_coef)
# train_env_scaled
# 
train_env_predictions_gas <- data.frame(predict(train_lm_gas, newdata = data.frame(train_ex_intermediate_gas)))
colnames(train_env_predictions_gas) <- c("predictions")
train_env_predictions_gas

lasso_dropcols <- names(elastic_coef_gas_terms[elastic_coef_gas_terms == 0])
env_clustering_terms_gas <- train_env_scaled_gas[,!(colnames(train_env_scaled_gas) %in% lasso_dropcols)]

env_clustering_terms_gas_pca <- prcomp(env_clustering_terms_gas, scale=FALSE)
gas_var <- get_pca_var(env_clustering_terms_gas_pca)

env_clustering_terms_gas_pca

tiff(
  "/Users/thomas/Desktop/Work/Research/uil/postquals/thermal/data/gas_pca.tiff",
  units="cm", 
  width=15, 
  height=15, 
  res=300, 
  pointsize=6
)

fviz_pca_var(
  env_clustering_terms_gas_pca,
  axes = c(1, 2),
  repel = TRUE,
  alpha.var = "contrib"
)

dev.off()

tiff(
  "/Users/thomas/Desktop/Work/Research/uil/postquals/thermal/data/gas_pca_d34.tiff",
  units="cm", 
  width=15, 
  height=15, 
  res=300, 
  pointsize=6
)

fviz_pca_var(
  env_clustering_terms_gas_pca,
  axes = c(3, 4),
  repel = TRUE,
  alpha.var = "contrib"
)

dev.off()

colnames(env_clustering_terms_gas) <- paste(colnames(env_clustering_terms_gas),"reg",sep="_")
exp(env_clustering_terms_gas) - 1

env_clustering_terms_gas

# we want to now just strip off the environmental variables, which this should do
env_terms_only_gas <- env_clustering_terms_gas[,2:(ncol(env_clustering_terms_gas)-3)]
env_terms_only_gas


train_envonly_predictions_gas <- data.frame(rowSums(env_terms_only_gas))
colnames(train_envonly_predictions_gas) <- c("env_predictions")

group_3_gas <- kmeans(env_terms_only_gas, centers=3)
group_5_gas <- kmeans(env_terms_only_gas, centers=5)
group_10_gas <- kmeans(env_terms_only_gas, centers=10)
# group_15 <- kmeans(env_terms_only, centers=15)
# group_20 <- kmeans(env_terms_only, centers=20)
# 
gmm_groups_gas$parameters$mean
abs(gmm_groups_gas$parameters$mean[,1] - gmm_groups_gas$parameters$mean[,2])

exp(sum(data.frame(gmm_groups_gas$parameters[2])[,"mean.1"])) - 1
exp(sum(data.frame(gmm_groups_gas$parameters[2])[,"mean.2"])) - 1

gmm_groups_gas <- Mclust(env_terms_only_gas, G=5)
gmm_groups_gas$parameters$mean
heatmap(gmm_groups_gas$parameters$mean)
# summary(gmm_groups_gas)
# plot(gmm_groups)
jpeg(
  "/Users/thomas/Desktop/Work/Research/uil/postquals/thermal/data/gas_gmm.jpeg",
  height=700,
  width=700,
  quality=300
)
# fviz_cluster(
#   gmm_groups_electric,
#   axes=c(2,3),
#   data=env_terms_only,
#   geom="point",
#   ggtheme = theme_light(),
#   pointsize=0.4,
#   dpi=300
# )

fviz_cluster(
  gmm_groups_gas,
  axes=c(1,2),
  data=env_terms_only_gas,
  geom="point",
  ggtheme = theme_light(),
  pointsize=0.4,
  dpi=300
)
dev.off()

# fviz_pca_ind(
#   env_terms_gas_pca, 
#   geom="point",
#   pointsize = 1,
#   col.ind = as.factor(group_3_gas$cluster),
#   palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#   addEllipses = TRUE, # Concentration ellipses
#   legend.title = "Groups",
#   select.ind = list(random = 500),
#   alpha=0.2
# )

# env_terms_gas_pca$rotation
# ggbiplot(
#   env_terms_gas_pca,
#   var.scale = 1,
#   groups = group_5_gas$cluster,
#   alpha = 0.05,
#   varname.size = 2
# )
# 
# 
# group_3_gas$centers

g5_cluster_gas <- as.matrix(gmm_groups_gas$classification)
colnames(g5_cluster_gas) <- c("cluster")

# predict(train_lm_gas, newdata = data.frame(train_ex_gas_t))
# # avg_profile <- colMeans(env_terms_only_gas)

# this is the estimated influence on the percentage of energy
100 * (exp(rowSums(group_5_gas$centers)) - 1)

percent_deviation_gas <- 100 * (exp(train_envonly_predictions_gas) - 1)
colnames(percent_deviation_gas) <- c("percent_deviation")

train_returning_km_gas <- cbind(
  gas_idx, 
  g5_cluster_gas, 
  (exp(env_terms_only_gas) - 1) * 100, 
  train_env_predictions_gas, 
  train_envonly_predictions_gas, 
  percent_deviation_gas
)

train_returning_km_gas
write.csv(train_returning_km_gas, "3c/training_environmental_gas_km5.csv", row.names = FALSE)

# gmm results
gmm_cluster_gas <- as.matrix(gmm_groups_gas$classification)
colnames(gmm_cluster_gas) <- c("cluster")

gmm_cluster_gas

train_returning_gmm_gas <- cbind(
  gas_idx, 
  gmm_cluster_gas, 
  env_terms_only_gas,
  train_env_predictions_gas, 
  train_envonly_predictions_gas, 
  percent_deviation_gas
)

write.csv(train_returning_gmm_gas, "3c/training_environmental_gas_gmm.csv", row.names = FALSE)
# 
# #### now want to try and run some simple regression models with what we have
# train_electric_pruned <- train_electric[,colnames(train_electric) %in% names(env_coef[env_coef != 0])]
# 
# train_electric_pruned$electric_norm
# #####
# 
# fit <- glm(electric_norm~., data=train_electric_pruned, family=gaussian())
# summary(fit)
# 
# #####
# train_gas <- data.frame(training_prepped) 
# train_gas$gas_norm <- gas_norm
# train_gas <- train_gas[!(train_gas$gas_norm == -Inf), ]
# 
# gas_fit <- glm(gas_norm~., data=train_gas, family=gaussian())
# summary(gas_fit)
# 
# 
# 
# 
# 
# 
# set.seed(42)
# cv_5 = trainControl(method = "cv", number = 5)
# 
# elnet = train(
#   electricity_mwh ~ ., train_electric = training_numeric_exits,
#   method = "glmnet",
#   trControl = cv_5
# )
# 
# x <- training_numeric_exits[, !names(training_numeric_exits) %in% c("coun_dist", "index_right","X1","bbl","bin","year","gas_mwh","electricity_mwh")]
# y <- training_numeric_exits$electricity_mwh
# 
# x_norm <- x %>% mutate_all(scale)
# 
# gm <- glmnet(x_norm, y)
# plot(gm)
# 
# cvfit <- cv.glmnet(as.matrix(x), y)
# plot(cvfit)
# 
# tmp_coeffs <- coef(cvfit, s = "lambda.min")
# coef_df <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
# nrow(coef_df)
# 
# coef_df[order(coef_df$coefficient),]
# 
# selected_names <- coef_df$name
# 
# glmnet_selected <- training_numeric_exits[, names(training_numeric_exits) %in% c(selected_names)]
# glmnet_selected$electricity_mwh <- y
# 
# simple_linear_model <- glm(electricity_mwh ~ ., family="gaussian", data=training_numeric)
# summary(simple_linear_model)
# 
# linear_model <- glm(electricity_mwh ~ ., family = "gaussian", data = glmnet_selected)
# summary(linear_model)


