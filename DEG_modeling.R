#DEGs_modeling.R
#This file takes in GSE174367 expression data and DEG file from GSE146639
#and applies ML models to classify AD vs. CTR

library(caret)
library(randomForest)
library(e1071)
library(readr)
library(pROC)
library(ggplot2)

#Replace with your project path to set base path for project
base_path <- "C:/Users/milla/OneDrive/Desktop/194b"
setwd(base_path)

#Load normalized GSE174367 expression data from GEO download
load("GSE174367_bulkRNA_processed.rda")
GSE174367_expr <- normExpr.reg

#Read in top 250 DEGs from GSE146639 microglial expression data
deg_features <- read_csv(file.path(base_path, "deg_features.csv"))

#Read in GSE174367 metadata and create AD/CTR labels
GSE174367_metadata <- read_excel(file.path(base_path, "GSE174367_metadata.xlsx"))
GSE174367_labels <- GSE174367_metadata$donor_group

#for reproducibility:
set.seed(42)

#Exploratory only- to inform later model decision
#Loop through different numbers of DEGS, classification model types
#choose counts and models to test:
deg_counts <- c(50, 100, 150, 200, 250)
models_to_test <- c("rf", "svmRadial", "glm")
results <- data.frame(NumGenes = integer(),
                      Model = character(),
                      Accuracy = numeric(),
                      stringsAsFactors = FALSE)
for (n in deg_counts) {
  # Subset top N DEGs
  top_n_genes <- deg_features$x[1:n]
  common_genes <- intersect(rownames(GSE174367_expr), top_n_genes)
  expr_subset <- GSE174367_expr[common_genes, , drop = FALSE]
  training_data <- as.data.frame(t(expr_subset))
  rownames(training_data) <- colnames(expr_subset)
  
  labels <- GSE174367_metadata$donor_group[match(rownames(training_data), GSE174367_metadata$Sample_id)]
  
  for (model in models_to_test) {
    if (model == "rf") {
      tune_grid <- expand.grid(mtry = c(2, 5, 10))
    } else if (model == "svmRadial") {
      tune_grid <- expand.grid(sigma = 0.01, C = c(0.25, 0.5, 1))
    } else {
      tune_grid <- NULL  
    }
    
    ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = multiClassSummary)
    
    tryCatch({
      trained_model <- train(
        x = training_data,
        y = as.factor(labels),
        method = model,
        trControl = ctrl,
        tuneGrid = tune_grid
      )
      best_acc <- max(trained_model$results$Accuracy)
      results <- add_row(results, NumGenes = n, Model = model, Accuracy = best_acc)
    }, error = function(e) {
      message(sprintf("Model '%s' with %d genes failed: %s", model, n, e$message))
    })
  }
}

print(results %>% arrange(desc(Accuracy)))

#Want to focus on RF-> performed well overall, useful for feature importance
#Baseline model: select top 150 DEGs
top_150_genes <- deg_features$x[1:150]
common_genes <- intersect(rownames(GSE174367_expr), top_150_genes)
expr_subset <- GSE174367_expr[common_genes, , drop = FALSE]

# Transpose: samples x genes
data_rf <- as.data.frame(t(expr_subset))
rownames(data_rf) <- colnames(expr_subset)

# Add labels (AD vs CTR)
labels_rf <- GSE174367_metadata$donor_group[match(rownames(data_rf), GSE174367_metadata$Sample_id)]
data_rf$Group <- as.factor(labels_rf)

#test train split
set.seed(42)
split_rf <- sample.split(data_rf$Group, SplitRatio = 0.7)
train_rf <- subset(data_rf, split_rf == TRUE)
test_rf <- subset(data_rf, split_rf == FALSE)

#train RF
set.seed(120)
classifier_RF <- randomForest(
  x = train_rf[, -ncol(train_rf)],  
  y = train_rf$Group,
  ntree = 500,
  importance = TRUE
)

#Evaluate RF
y_pred_class <- predict(classifier_RF, newdata = test_rf[, -ncol(test_rf)])

# Create confusion matrix
conf_mat <- confusionMatrix(y_pred_class, test_rf$Group)
accuracy_rf <- conf_mat$overall["Accuracy"]
precision_rf <- conf_mat$byClass["Precision"]
recall_rf <- conf_mat$byClass["Recall"]

# Predict probabilities
y_prob <- predict(classifier_RF, newdata = test_rf[, -ncol(test_rf)], type = "prob")[, 2]

# ROC and AUC
roc_obj <- roc(response = test_rf$Group, predictor = y_prob, levels = rev(levels(test_rf$Group)))
auroc_rf <- auc(roc_obj)

cat("Accuracy :", round(accuracy_rf, 4), "\n")
cat("Precision:", round(precision_rf, 4), "\n")
cat("Recall   :", round(recall_rf, 4), "\n")
cat("AUROC    :", round(auroc_rf, 4), "\n")

#feature importance to get top 10 important genes
importance_scores <- importance(classifier_RF)
importance_df <- data.frame(Gene = rownames(importance_scores), Importance = importance_scores[, 1])
importance_df <- importance_df[order(-importance_df$Importance), ]
print(head(importance_df, 10))


# Plot the top 10 important genes using ggplot2
library(ggplot2)
ggplot(top_10_importance, aes(x = reorder(Gene, Importance), y = Importance)) +
  geom_point(color = "steelblue", size = 4) +  
  coord_flip() +  
  labs(title = "Top 10 Important Genes (RF)", x = "Gene", y = "Importance") +
  theme_minimal()

#Plot ROC curve
plot(roc_obj, main = "ROC Curve - Final RF Model", col = "blue", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray")
