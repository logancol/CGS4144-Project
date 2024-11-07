#install.packages("cluster")
library(cluster)

# library init
library(data.table)
library(tidyclust)
#install.packages("tidymodels")
library(tidymodels)
#install.packages("mclust")
library(mclust)
library(dplyr)
# df init
data <- fread("HUGO_Final.tsv")
differentially_expressed <- fread("differentially_expressed_genes.csv")
metadata <- fread("metadata_SRP192714.tsv")

# getting top 5000 highest variance genes

# saving hugo names (I'm going to cluster first using Ensembl)
hugo_names <- data %>%
  select(HUGO, Ensembl)
rownames(data) <- data$Ensembl
data <- data %>% select (-c( HUGO))
expression_data <- data[, -1, with=FALSE]
variances <- apply(expression_data, 1, var)
sorted_variances <- sort(variances, decreasing = TRUE)
top_5000_indices <- order(variances, decreasing = TRUE)[1:5000]
top_5000_genes <- rownames(data)[top_5000_indices]
top_5000_data <- data[top_5000_indices, ]
top_5000_data <- top_5000_data[,-1]
top_5000_data <- t(top_5000_data)
top_10_data <- top_5000_data[,1:10]
top_100_data <- top_5000_data[,1:100]
top_1000_data <- top_5000_data[,1:1000]
genes_10000 <- t(expression_data)
genes_10000 <- genes_10000[,1:10000]

outcome <- data.frame(metadata$refinebio_accession_code,metadata$refinebio_time)
rownames(outcome) <- outcome$metadata.refinebio_accession_code 
outcome <- outcome[,-1]

install.packages("glmnet")
library(glmnet)

# Prepare data for glmnet (excluding the outcome column)
x_5000 <- as.matrix(top_5000_data)
y_5000 <- as.factor(outcome)

# Fit a regularized multinomial logistic regression model with glmnet
cv_fit1 <- cv.glmnet(x_5000, y_5000, family = "multinomial", alpha = 0.5)  # alpha=0.5 for elastic-net
plot(cv_fit1)

#clustering
pam_result <- pam(top_5000_data, k = 3)
cluster_assignments <- pam_result$clustering
x_cluster <- as.matrix(top_5000_data)  # Predictor matrix
y_cluster <- cluster_assignments  # Outcome (target) variable
cv_fit2 <- cv.glmnet(x_cluster, y_cluster, family = "multinomial", alpha = 0.5)
plot(cv_fit2)

#10 genes
x_10 <- as.matrix(top_10_data)
y_10 <- as.factor(outcome)
cv_fit3 <- cv.glmnet(x_10, y_10, family = "multinomial", alpha = 0.5)  # alpha=0.5 for elastic-net
plot(cv_fit3)

#100 genes
x_100 <- as.matrix(top_100_data)
y_100 <- as.factor(outcome)
cv_fit4 <- cv.glmnet(x_100, y_100, family = "multinomial", alpha = 0.5)  # alpha=0.5 for elastic-net
plot(cv_fit4)

#1000 genes
x_1000 <- as.matrix(top_1000_data)
y_1000 <- as.factor(outcome)
cv_fit5 <- cv.glmnet(x_1000, y_1000, family = "multinomial", alpha = 0.5)  # alpha=0.5 for elastic-net
plot(cv_fit5)

#10000 genes
x_10000 <- as.matrix(genes_10000)
y_10000 <- as.factor(outcome)
cv_fit6 <- cv.glmnet(x_10000, y_10000, family = "multinomial", alpha = 0.5)  # alpha=0.5 for elastic-net
plot(cv_fit6)

#data1 <- data.frame(metadata$refinebio_accession_code,cluster_assignments,metadata$refinebio_time)

#confusion matrix function
generate_confusion_matrix <- function(model, predictors, true_labels) {
  predictions <- predict(model, newx = predictors, s = "lambda.min", type = "class")
  confusion_matrix <- table(True = true_labels, Predicted = predictions)
  return(confusion_matrix)
}

confusion_matrix_5000 <- generate_confusion_matrix(cv_fit1, x_5000, y_5000)
print(confusion_matrix_5000)


# Install and load the pROC package
#install.packages("pROC")
library(pROC)

# Get predicted probabilities for each class
predicted_probs_5000 <- predict(cv_fit1, newx = x_5000, s = "lambda.min", type = "response")
predicted_probs_5000 <- data.frame(predicted_probs_5000)
colnames(predicted_probs_5000) <- c("convalescent", "early acute", "late acute")
y_5000 <- factor(y_5000)
# Calculate AUC for each class vs. rest using multiclass.roc
multiclass_roc <- multiclass.roc(y_5000, predicted_probs_5000)

predicted_probs_10 <- predict(cv_fit3, newx = x_10, s = "lambda.min", type = "response")
predicted_probs_10 <- data.frame(predicted_probs_10)
colnames(predicted_probs_10) <- c("convalescent", "early acute", "late acute")
y_10 <- factor(y_10)
multiclass_roc10 <- multiclass.roc(y_10, predicted_probs_10)

predicted_probs_100 <- predict(cv_fit4, newx = x_100, s = "lambda.min", type = "response")
predicted_probs_100 <- data.frame(predicted_probs_100)
colnames(predicted_probs_100) <- c("convalescent", "early acute", "late acute")
y_100 <- factor(y_100)
multiclass_roc100 <- multiclass.roc(y_100, predicted_probs_100)

predicted_probs_1000 <- predict(cv_fit5, newx = x_1000, s = "lambda.min", type = "response")
predicted_probs_1000 <- data.frame(predicted_probs_1000)
colnames(predicted_probs_1000) <- c("convalescent", "early acute", "late acute")
y_1000 <- factor(y_1000)
multiclass_roc1000 <- multiclass.roc(y_1000, predicted_probs_1000)

predicted_probs_10000 <- predict(cv_fit6, newx = x_10000, s = "lambda.min", type = "response")
predicted_probs_10000 <- data.frame(predicted_probs_10000)
colnames(predicted_probs_10000) <- c("convalescent", "early acute", "late acute")
y_10000 <- factor(y_10000)
multiclass_roc10000 <- multiclass.roc(y_10000, predicted_probs_10000)



# Print AUC result

cat(" 10 genes AUC:", multiclass_roc10$auc, "\n",
    "100 genes AUC:", multiclass_roc100$auc, "\n", "1000 genes AUC:", multiclass_roc1000$auc, "\n",
    "10000 genes AUC:", multiclass_roc10000$auc, "\n")

# accuracy
calculate_accuracy <- function(model, predictors, true_labels) {
  predicted_classes <- predict(model, newx = predictors, s = "lambda.min", type = "class")
  accuracy <- mean(predicted_classes == true_labels)
  return(accuracy)
}

# Calculate accuracy for each model
accuracy_10 <- calculate_accuracy(cv_fit3, x_10, y_10)
accuracy_100 <- calculate_accuracy(cv_fit4, x_100, y_100)
accuracy_1000 <- calculate_accuracy(cv_fit5, x_1000, y_1000)
accuracy_5000 <- calculate_accuracy(cv_fit1, x_5000, y_5000)
accuracy_10000 <- calculate_accuracy(cv_fit6, x_10000, y_10000)

# Print the accuracies
cat("Accuracy for top 5000 genes:", accuracy_5000, "\n", "5000 genes AUC:", multiclass_roc$auc, "\n")

cat(" Accuracy for top 10 genes:", accuracy_10,"\n", "Accuracy for top 100 genes:", accuracy_100,"\n",
    "Accuracy for top 1000 genes:", accuracy_1000,"\n", "Accuracy for top 10000 genes:", accuracy_10000, "\n")




