# Install and load required packages
install.packages(c("tidymodels", "ranger", "ggplot2", "reshape2", "data.table"))
library(tidymodels)
library(ranger)
library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)
library(caret)

# Load data
expression_data <- fread("HUGO_Final.tsv")  # Expression data
cluster_data <- read.csv("consensus.csv", row.names = 1)  # Cluster assignments
metadata <- read.delim("metadata_SRP192714.tsv")  # Metadata

# Print initial diagnostics
print("Initial data checks:")
print("Expression data dimensions:")
print(dim(expression_data))
print("\nCluster data dimensions:")
print(dim(cluster_data))
print("\nNumber of clusters:")
print(table(cluster_data$x))

# Convert expression data to numeric matrix
gene_names <- expression_data$Ensembl
expression_matrix <- as.matrix(expression_data[, -c("Ensembl", "HUGO"), with=FALSE])
rownames(expression_matrix) <- gene_names

# Verify data is numeric
print("\nChecking data types:")
print("Is matrix numeric?")
print(is.numeric(expression_matrix))
print("\nAny NA values?")
print(sum(is.na(expression_matrix)))

# Calculate top 5000 most variable genes
variances <- apply(expression_matrix, 1, var)
sorted_variances <- sort(variances, decreasing = TRUE)
top_5000_indices <- order(variances, decreasing = TRUE)[1:5000]
top_5000_data <- expression_matrix[top_5000_indices, ]

# Transpose and scale the data
expression_data_t <- t(top_5000_data)
expression_data_t <- scale(expression_data_t)
expression_data_t <- as.data.frame(expression_data_t)

# Match cluster assignments
cluster_assignments <- cluster_data[rownames(expression_data_t), "x"]
print("\nNumber of NA values in cluster assignments:")
print(sum(is.na(cluster_assignments)))

# Add cluster labels
expression_data_t$Label <- factor(cluster_assignments)

# Set up Random Forest model
rf_spec <- 
  rand_forest(
    mtry = floor(sqrt(ncol(expression_data_t) - 1)),
    trees = 500,
    min_n = 2
  ) %>%
  set_mode("classification") %>%
  set_engine("ranger", importance = "impurity")

# Split data into training and testing sets
set.seed(123)
data_split <- initial_split(expression_data_t, prop = 0.7)
train_data <- training(data_split)
test_data <- testing(data_split)

# Train initial model with 5000 genes
rf_fit <- rf_spec %>% 
  fit(Label ~ ., data = train_data)

# Make predictions
predictions <- predict(rf_fit, test_data)
results <- tibble(
  truth = test_data$Label,
  prediction = predictions$.pred_class
)

# Make predictions for all samples
all_rf_predictions_cluster <- predict(rf_fit, expression_data_t)

# Create a new RF model for predicting infection stages
expression_data_t_stages <- expression_data_t
expression_data_t_stages$Label <- metadata$refinebio_time[match(rownames(expression_data_t), 
                                                                metadata$refinebio_accession_code)]
expression_data_t_stages$Label <- factor(expression_data_t_stages$Label)

# Train new RF model for stages
set.seed(123)
rf_fit_stages <- rf_spec %>% 
  fit(Label ~ ., data = expression_data_t_stages)

# Make predictions for infection stages
all_rf_predictions_stages <- predict(rf_fit_stages, expression_data_t_stages)

# Format predictions for all samples
rf_formatted_all <- data.frame(
  Sample = rownames(expression_data_t),
  Cluster_5000_RF = cluster_assignments,
  Model_5000_RF = as.character(all_rf_predictions_cluster$.pred_class),
  Time_5000_RF = as.character(all_rf_predictions_stages$.pred_class)
)

# Save complete RF predictions
write.csv(rf_formatted_all, "rf_predictions_all_samples.csv", row.names = FALSE)

# Get feature importance
importance_scores <- rf_fit$fit$variable.importance
importance_df <- data.frame(
  feature = names(importance_scores),
  importance = importance_scores
) %>%
  arrange(desc(importance))

# Save top 20 important features
write.csv(head(importance_df, 20), "top_20_important_features.csv")

# Test different numbers of genes
gene_counts <- c(10, 100, 1000, 5000)
accuracies <- data.frame(genes = gene_counts, accuracy = NA)

print("\nTesting different numbers of genes:")
for(i in seq_along(gene_counts)) {
  n_genes <- gene_counts[i]
  print(paste("Processing", n_genes, "genes..."))
  
  # Get top n genes
  current_indices <- order(variances, decreasing = TRUE)[1:n_genes]
  current_data <- expression_matrix[current_indices, ]
  
  # Transpose and scale
  current_data_t <- t(current_data)
  current_data_t <- scale(current_data_t)
  current_data_t <- as.data.frame(current_data_t)
  current_data_t$Label <- factor(cluster_assignments)
  
  # Update RF spec for current feature set
  current_rf_spec <- rand_forest(
    mtry = floor(sqrt(n_genes)),
    trees = 500,
    min_n = 2
  ) %>%
    set_mode("classification") %>%
    set_engine("ranger", importance = "impurity")
  
  # Split data
  set.seed(123)
  current_split <- initial_split(current_data_t, prop = 0.7)
  current_train <- training(current_split)
  current_test <- testing(current_split)
  
  # Fit model and predict
  current_fit <- current_rf_spec %>% fit(Label ~ ., data = current_train)
  current_pred <- predict(current_fit, current_test)
  
  # Calculate accuracy
  accuracies$accuracy[i] <- mean(current_pred$.pred_class == current_test$Label)
  
  # Get and save feature importance
  current_importance <- data.frame(
    feature = names(current_fit$fit$variable.importance),
    importance = current_fit$fit$variable.importance
  ) %>%
    arrange(desc(importance))
  
  write.csv(current_importance, 
            paste0("rf_feature_importance_", n_genes, "_genes.csv"), 
            row.names = FALSE)
  
  # Print confusion matrix
  current_conf <- conf_mat(
    data = tibble(
      truth = current_test$Label,
      prediction = current_pred$.pred_class
    ),
    truth,
    prediction
  )
  print(paste("\nConfusion Matrix for", n_genes, "genes:"))
  print(current_conf)
}

# Create and save accuracy plot
acc_plot <- ggplot(accuracies, aes(x = genes, y = accuracy)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  labs(title = "Random Forest Model Accuracy vs Number of Genes",
       x = "Number of Genes",
       y = "Accuracy") +
  theme_minimal()

# Save results
write.csv(accuracies, "rf_cluster_prediction_results.csv")
ggsave("rf_cluster_prediction_accuracy.png", acc_plot)

# Binary classification (early acute vs others)
gene_counts <- c(10, 100, 1000, 10000)
results_df <- data.frame(
  genes = gene_counts,
  accuracy = numeric(length(gene_counts)),
  auc = numeric(length(gene_counts))
)

for(i in seq_along(gene_counts)) {
  n_genes <- gene_counts[i]
  print(paste("Processing", n_genes, "genes for binary classification..."))
  
  current_indices <- order(variances, decreasing = TRUE)[1:n_genes]
  current_data <- expression_matrix[current_indices, ]
  
  current_data_t <- t(current_data)
  current_data_t <- scale(current_data_t)
  current_data_t <- as.data.frame(current_data_t)
  
  # Create binary labels (early acute vs others)
  binary_labels <- factor(ifelse(
    metadata$refinebio_time[match(rownames(current_data_t), 
                                  metadata$refinebio_accession_code)] == "early acute",
    "early_acute", "others"
  ))
  
  current_data_t$Label <- binary_labels
  
  # Update RF spec for current feature set
  current_rf_spec <- rand_forest(
    mtry = floor(sqrt(n_genes)),
    trees = 500,
    min_n = 2
  ) %>%
    set_mode("classification") %>%
    set_engine("ranger", importance = "impurity")
  
  set.seed(123)
  current_split <- initial_split(current_data_t, prop = 0.7)
  current_train <- training(current_split)
  current_test <- testing(current_split)
  
  current_fit <- current_rf_spec %>% fit(Label ~ ., data = current_train)
  current_pred <- predict(current_fit, current_test)
  current_pred_prob <- predict(current_fit, current_test, type = "prob")
  
  results_df$accuracy[i] <- mean(current_pred$.pred_class == current_test$Label)
  
  results_df$auc[i] <- roc_auc(data.frame(
    truth = current_test$Label,
    pred = current_pred_prob$.pred_early_acute
  ), truth, pred)$.estimate
  
  # Save feature importance
  current_importance <- data.frame(
    feature = names(current_fit$fit$variable.importance),
    importance = current_fit$fit$variable.importance
  ) %>%
    arrange(desc(importance))
  write.csv(current_importance, 
            paste0("rf_binary_feature_importance_", n_genes, "_genes.csv"), 
            row.names = FALSE)
}

# Create performance plot
performance_plot <- ggplot(results_df %>% 
                             gather(key = "metric", value = "value", -genes)) +
  geom_line(aes(x = genes, y = value, color = metric)) +
  geom_point(aes(x = genes, y = value, color = metric)) +
  scale_x_log10() +
  labs(title = "Random Forest Performance vs Number of Genes",
       x = "Number of Genes",
       y = "Value",
       color = "Metric") +
  theme_minimal()

# Save results
ggsave("rf_performance_vs_genes.png", performance_plot)
write.csv(results_df, "rf_binary_classification_results.csv", row.names = FALSE)

# Create and save summary
summary_text <- paste(
  "Random Forest Analysis Results:\n",
  "\nMain Classification Performance:",
  paste("Accuracy with 5000 genes:", 
        round(mean(predictions$.pred_class == test_data$Label), 3)),
  "\nBinary Classification (Early Acute vs Others) Results:",
  "\nPerformance by gene count:",
  paste("10 genes - Accuracy:", round(results_df$accuracy[1], 3), 
        "AUC:", round(results_df$auc[1], 3)),
  paste("100 genes - Accuracy:", round(results_df$accuracy[2], 3), 
        "AUC:", round(results_df$auc[2], 3)),
  paste("1000 genes - Accuracy:", round(results_df$accuracy[3], 3), 
        "AUC:", round(results_df$auc[3], 3)),
  paste("10000 genes - Accuracy:", round(results_df$accuracy[4], 3), 
        "AUC:", round(results_df$auc[4], 3)),
  "\nFeature importance files have been saved for each analysis.",
  sep = "\n"
)

# Print and save summary
cat(summary_text)
writeLines(summary_text, "rf_analysis_summary.txt")

# Print final summary
print("\nAnalysis complete. Files saved:")
print("1. rf_predictions_all_samples.csv - All predictions")
print("2. top_20_important_features.csv - Top important features")
print("3. rf_cluster_prediction_results.csv - Accuracy results")
print("4. rf_binary_classification_results.csv - Binary classification results")
print("5. Feature importance files for each gene count")
print("6. Visualization plots")


# Load required libraries
library(tidymodels)
library(ggplot2)
library(caret)
library(viridis)

# Create confusion matrices for both cluster and time predictions
# For cluster predictions
cluster_conf <- confusionMatrix(
  data = factor(rf_formatted_all$Model_5000_RF),
  reference = factor(rf_formatted_all$Cluster_5000_RF)
)

# For time predictions
time_conf <- confusionMatrix(
  data = factor(rf_formatted_all$Time_5000_RF),
  reference = factor(metadata$refinebio_time[match(rf_formatted_all$Sample, 
                                                   metadata$refinebio_accession_code)])
)

# Function to create confusion matrix plot
plot_confusion_matrix <- function(conf_matrix, title) {
  # Convert to data frame
  conf_data <- as.data.frame(conf_matrix$table)
  
  # Calculate percentages
  conf_data$Percentage <- conf_data$Freq / sum(conf_data$Freq) * 100
  
  # Create plot
  ggplot(conf_data, aes(x = Reference, y = Prediction, fill = Percentage)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Percentage, Freq)), 
              color = "black", size = 3) +
    scale_fill_viridis() +
    labs(title = title,
         x = "True Label",
         y = "Predicted Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Create plots
cluster_plot <- plot_confusion_matrix(cluster_conf, "Confusion Matrix - Cluster Predictions")
time_plot <- plot_confusion_matrix(time_conf, "Confusion Matrix - Time Point Predictions")

# Save plots
ggsave("cluster_confusion_matrix.png", cluster_plot, width = 8, height = 6)
ggsave("time_confusion_matrix.png", time_plot, width = 10, height = 8)

# Print detailed statistics
print("Cluster Prediction Statistics:")
print(cluster_conf)
print("\nTime Point Prediction Statistics:")
print(time_conf)

# Create summary statistics
cluster_summary <- data.frame(
  Metric = c("Accuracy", "Kappa", "Balanced Accuracy"),
  Value = c(
    cluster_conf$overall["Accuracy"],
    cluster_conf$overall["Kappa"],
    mean(diag(prop.table(cluster_conf$table, 1)))
  )
)

time_summary <- data.frame(
  Metric = c("Accuracy", "Kappa", "Balanced Accuracy"),
  Value = c(
    time_conf$overall["Accuracy"],
    time_conf$overall["Kappa"],
    mean(diag(prop.table(time_conf$table, 1)))
  )
)

# Save summaries
write.csv(cluster_summary, "cluster_prediction_metrics.csv", row.names = FALSE)
write.csv(time_summary, "time_prediction_metrics.csv", row.names = FALSE)

# Calculate per-class metrics for clusters
cluster_per_class <- data.frame(
  Class = rownames(cluster_conf$byClass),
  Sensitivity = cluster_conf$byClass[, "Sensitivity"],
  Specificity = cluster_conf$byClass[, "Specificity"],
  Precision = cluster_conf$byClass[, "Pos Pred Value"],
  F1_Score = cluster_conf$byClass[, "F1"]
)

# Calculate per-class metrics for time points
time_per_class <- data.frame(
  Class = rownames(time_conf$byClass),
  Sensitivity = time_conf$byClass[, "Sensitivity"],
  Specificity = time_conf$byClass[, "Specificity"],
  Precision = time_conf$byClass[, "Pos Pred Value"],
  F1_Score = time_conf$byClass[, "F1"]
)

# Save per-class metrics
write.csv(cluster_per_class, "cluster_per_class_metrics.csv", row.names = FALSE)
write.csv(time_per_class, "time_per_class_metrics.csv", row.names = FALSE)

# Create summary text
summary_text <- paste(
  "Confusion Matrix Analysis Summary\n",
  "\nCluster Predictions:",
  paste("Overall Accuracy:", round(cluster_conf$overall["Accuracy"], 3)),
  paste("Kappa:", round(cluster_conf$overall["Kappa"], 3)),
  "\nPer-class Performance (Sensitivity):",
  paste(capture.output(round(cluster_conf$byClass[, "Sensitivity"], 3)), collapse = "\n"),
  "\n\nTime Point Predictions:",
  paste("Overall Accuracy:", round(time_conf$overall["Accuracy"], 3)),
  paste("Kappa:", round(time_conf$overall["Kappa"], 3)),
  "\nPer-class Performance (Sensitivity):",
  paste(capture.output(round(time_conf$byClass[, "Sensitivity"], 3)), collapse = "\n"),
  sep = "\n"
)

# Save summary
writeLines(summary_text, "confusion_matrix_summary.txt")