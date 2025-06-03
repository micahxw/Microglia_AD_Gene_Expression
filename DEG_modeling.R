#DEGs_modeling.R
#This file takes in GSE174367 expression data and DEG file from GSE146639
#and applies ML models to classify AD vs. CTR

library(readxl)
library(readr)
library(dplyr)
library(caret)
library(pROC)
library(org.Hs.eg.db)
library(ggplot2)

#Replace with your project path to set base path for project
base_path <- "path/to/project/folder"
setwd(base_path)

#Load normalized GSE174367 expression data from GEO download
load("GSE174367_bulkRNA_processed.rda")
GSE174367_expr <- normExpr.reg

#Load metadata
GSE174367_metadata <- read_excel(file.path(base_path, "GSE174367_metadata.xlsx"))

#Read in top DEGs from GSE146639 microglial expression data
deg_features <- read_csv(file.path(base_path, "deg_features.csv"))


#use top 60 DEGS, subset GSE174367 to these genes
top_genes <- deg_features$x[1:60]
common_genes <- intersect(rownames(GSE174367_expr), top_genes)
expr_subset <- GSE174367_expr[common_genes, , drop = FALSE]

#prep training data and labels
training_data <- as.data.frame(t(expr_subset))
rownames(training_data) <- colnames(expr_subset)
labels <- GSE174367_metadata$donor_group[match(rownames(training_data), GSE174367_metadata$Sample_id)]
labels <- factor(labels, levels = c("CTR", "AD"))

#evaluation function for repeated model metrics
evaluate_model <- function(pred, prob, actual){
  actual <- factor(actual, levels = c("CTR", "AD"))
  pred <- factor(pred, levels = c("CTR", "AD"))
  
  cm <- confusionMatrix(pred, actual, positive="AD")
  auc <- roc(actual, prob, levels = c("CTR", "AD"))$auc
  data.frame(
    Accuracy = cm$overall["Accuracy"],
    Precision = cm$byClass["Precision"],
    Recall = cm$byClass["Recall"],
    F1 = cm$byClass["F1"],
    AUROC = auc
  )
}
set.seed(123)

# set up 5-fold CV with 3 repeats
cv <- trainControl(method='repeatedcv', number=5, classProbs=TRUE, repeats=3,
                   summaryFunction = twoClassSummary,savePredictions = "final")

#RF(baseline) model
rf_model<- train(x = training_data, y=labels,
                 method="rf", 
                 metric="ROC", 
                 trControl=cv)


rf_pred <- rf_model$pred
rf_pred_labels <- factor(ifelse(rf_pred$AD >= 0.5, "AD", "CTR"), levels = c("CTR", "AD"))
rf_metrics <- evaluate_model(rf_pred_labels, rf_pred$AD, rf_pred$obs)
print(rf_metrics)


#XGBOOST (baseline) model
xgb_model <- train(x = training_data, y = labels, 
  method = "xgbTree", 
  metric = "ROC", 
  trControl = cv)


xgb_preds <- xgb_model$pred
xgb_pred_labels <- factor(ifelse(xgb_preds$AD >= 0.5, "AD", "CTR"), levels = c("CTR", "AD"))
xgb_metrics <- evaluate_model(xgb_pred_labels, xgb_preds$AD, xgb_preds$obs)
print(xgb_metrics)


#LASSO- tuning + feature importance
tune_grid <- expand.grid(alpha = 1,lambda = 10^seq(-4, 1, length = 20))

lasso_model <- train(x = training_data, y = labels,
                      method = "glmnet",
                      metric = "ROC",
                      trControl = cv,
                      tuneGrid = tune_grid,
                     preProcess = c("center", "scale"))


best_lambda <- lasso_model$bestTune$lambda
lasso_preds <- lasso_model$pred %>% filter(lambda == best_lambda)
roc_obj <- roc(lasso_preds$obs, lasso_preds$AD, levels = c("CTR", "AD"))

#find best threshold for LASSO metrics
best_threshold <- as.numeric(coords(roc_obj, "best", ret = "threshold"))
print(best_threshold)
lasso_pred_labels <- factor(ifelse(lasso_preds$AD >= best_threshold, "AD", "CTR"), levels = c("CTR", "AD"))
lasso_metrics <- evaluate_model(lasso_pred_labels, lasso_preds$AD, lasso_preds$obs)
print(lasso_metrics)

# Feature Importance from LASSO
# Extract coefficients for best lambda, arrange from best to worst
coef_best <- coef(lasso_model$finalModel, s = best_lambda)
feature_importance <- data.frame(
  Gene = rownames(coef_best),
  Coefficient = as.vector(coef_best))%>% 
  filter(Gene != "(Intercept)")%>%
  arrange(desc(abs(Coefficient)))


# Get top 10 features by absolute coefficient
top_features <- feature_importance %>%
  slice_max(order_by = abs(Coefficient), n = 10)

#convert features from ENSEMBL ID to gene names for visualization
map_gene_names <- mapIds(org.Hs.eg.db, keys=top_features$Gene, column='SYMBOL',
                         keytype = "ENSEMBL")
top_features$gene_names <- map_gene_names

# Plot in dot plot style
ggplot(top_features, aes(x = reorder(gene_names, abs(Coefficient)), y = Coefficient)) +
  geom_point(aes(color = Coefficient > 0), size = 4) +
  scale_color_manual(values = c("steelblue", "firebrick"), guide = "none") +
  coord_flip() +
  labs(title = "Top 10 Important Genes in AD vs. CTR Classification",
       x = "Gene",
       y = "Lasso Coefficient") 

