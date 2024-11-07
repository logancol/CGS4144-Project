# Install and load required packages
install.packages(c("kknn", "tidymodels", "ggplot2", "reshape2", "data.table"))
library(tidymodels)
library(kknn)
library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)

# Load data
# Note: Ensure your files are in the working directory
expression_data <- fread("HUGO_Final.tsv")  # Your expression data
cluster_data <- read.delim("for_heatmap.tsv", sep="\t")  # Your cluster assignments

# Print initial diagnostics
print("Initial data checks:")
print("Expression data dimensions:")
print(dim(expression_data))
print("\nNumber of clusters:")
print(table(cluster_data$Cluster_5000))

# Convert expression data to numeric matrix
# First, store gene names
gene_names <- expression_data$Ensembl
# Remove non-numeric columns and convert to matrix
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
expression_data_t <- scale(expression_data_t)  # Scale numeric data
expression_data_t <- as.data.frame(expression_data_t)  # Convert back to dataframe

# Match cluster assignments with expression data
cluster_assignments <- cluster_data$Cluster_5000[match(rownames(expression_data_t), 
                                                       cluster_data$Sample)]
expression_data_t$Label <- factor(cluster_assignments)

# Set up KNN model
knn_spec <- 
  nearest_neighbor(neighbors = 5) %>%
  set_mode("classification") %>%
  set_engine("kknn")

# Split data into training and testing sets
set.seed(123)  # For reproducibility
data_split <- initial_split(expression_data_t, prop = 0.7)
train_data <- training(data_split)
test_data <- testing(data_split)

# Train initial model with 5000 genes
knn_fit <- knn_spec %>% 
  fit(Label ~ ., data = train_data)

# Make predictions
predictions <- predict(knn_fit, test_data)
results <- tibble(
  truth = test_data$Label,
  prediction = predictions$.pred_class
)

# Make predictions for cluster assignments
all_knn_predictions_cluster <- predict(knn_fit, expression_data_t)

# Create a new KNN model for predicting infection stages
# First, create a new data frame with infection stages as labels
expression_data_t_stages <- expression_data_t
expression_data_t_stages$Label <- metadata$refinebio_time[match(rownames(expression_data_t), 
                                                                metadata$refinebio_accession_code)]
expression_data_t_stages$Label <- factor(expression_data_t_stages$Label)

# Train new KNN model for stages
set.seed(123)
knn_fit_stages <- knn_spec %>% 
  fit(Label ~ ., data = expression_data_t_stages)

# Make predictions for infection stages
all_knn_predictions_stages <- predict(knn_fit_stages, expression_data_t_stages)

# Format predictions for all samples
knn_formatted_all <- data.frame(
  Sample = rownames(expression_data_t),
  Cluster_5000_KNN = cluster_assignments,
  Model_5000_KNN = as.character(all_knn_predictions_cluster$.pred_class),
  Time_5000_KNN = as.character(all_knn_predictions_stages$.pred_class)
)

# Save complete KNN predictions
write.csv(knn_formatted_all, "knn_predictions_all_samples.csv", row.names = FALSE)

metadata <- read.delim("metadata_SRP192714.tsv")
stage_mapping <- metadata$refinebio_time[match(cluster_data$Sample, metadata$refinebio_accession_code)]
names(stage_mapping) <- cluster_data$Cluster_5000

# Now create the formatted predictions
knn_formatted <- data.frame(
  Sample = rownames(test_data),
  Cluster_5000_KNN = test_data$Label,
  Model_5000_KNN = stage_mapping[as.character(predictions$.pred_class)]  # Convert predictions to stage names
)

write.csv(knn_formatted, "knn_predictions_formatted.csv", row.names = FALSE)


# Generate and print confusion matrix
print("\nConfusion Matrix for 5000 genes:")
conf_matrix <- conf_mat(results, truth, prediction)
print(conf_matrix)

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
  
  # Add labels
  current_data_t$Label <- factor(cluster_assignments)
  
  # Split data
  set.seed(123)
  current_split <- initial_split(current_data_t, prop = 0.7)
  current_train <- training(current_split)
  current_test <- testing(current_split)
  
  # Fit model and predict
  current_fit <- knn_spec %>% fit(Label ~ ., data = current_train)
  current_pred <- predict(current_fit, current_test)
  
  # Calculate accuracy
  current_accuracy <- mean(current_pred$.pred_class == current_test$Label)
  accuracies$accuracy[i] <- current_accuracy
  
  # Print results for this iteration
  print(paste("Accuracy with", n_genes, "genes:", round(current_accuracy, 3)))
  
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
  labs(title = "KNN Model Accuracy vs Number of Genes",
       x = "Number of Genes",
       y = "Accuracy") +
  theme_minimal()

# Save results
print("\nSaving results...")
write.csv(accuracies, "knn_cluster_prediction_results.csv")
write.csv(as.data.frame(conf_matrix$table), "knn_confusion_matrix.csv")
ggsave("cluster_prediction_accuracy.png", acc_plot)

# Print final summary
print("\nFinal Results Summary:")
print("1. Accuracy with different numbers of genes:")
print(accuracies)
print("\n2. Confusion Matrix for final model:")
print(conf_matrix)
print("\n3. Cluster distribution:")
print(table(cluster_assignments))


library(caret)

# For KNN:
# Confusion matrix for cluster predictions
knn_cluster_conf <- confusionMatrix(
  data = factor(knn_formatted_all$Model_5000_KNN),
  reference = factor(knn_formatted_all$Cluster_5000_KNN)
)
print("KNN Cluster Prediction Confusion Matrix:")
print(knn_cluster_conf)

# Confusion matrix for stage predictions
knn_stage_conf <- confusionMatrix(
  data = factor(knn_formatted_all$Time_5000_KNN),
  reference = factor(metadata$refinebio_time[match(knn_formatted_all$Sample, metadata$refinebio_accession_code)])
)
print("\nKNN Stage Prediction Confusion Matrix:")
print(knn_stage_conf)




# Set up gene counts to test
gene_counts <- c(10, 100, 1000, 10000)
results_df <- data.frame(
  genes = gene_counts,
  accuracy = numeric(length(gene_counts)),
  auc = numeric(length(gene_counts))  # Single AUC for binary classification
)

# For each gene count
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
  
  # Create binary labels (early acute vs others)
  binary_labels <- factor(ifelse(
    metadata$refinebio_time[match(rownames(current_data_t), metadata$refinebio_accession_code)] == "early acute",
    "early_acute", "others"
  ))
  
  current_data_t$Label <- binary_labels
  
  # Split data
  set.seed(123)
  current_split <- initial_split(current_data_t, prop = 0.7)
  current_train <- training(current_split)
  current_test <- testing(current_split)
  
  # Fit model and predict
  current_fit <- knn_spec %>% fit(Label ~ ., data = current_train)
  current_pred <- predict(current_fit, current_test)
  current_pred_prob <- predict(current_fit, current_test, type = "prob")
  
  # Calculate accuracy
  results_df$accuracy[i] <- mean(current_pred$.pred_class == current_test$Label)
  
  # Calculate AUC for binary classification
  results_df$auc[i] <- roc_auc(data.frame(
    truth = current_test$Label,
    pred = current_pred_prob$.pred_early_acute
  ), truth, pred)$.estimate
  
  # Print confusion matrix
  conf_matrix <- conf_mat(
    data = tibble(
      truth = current_test$Label,
      prediction = current_pred$.pred_class
    ),
    truth,
    prediction
  )
  print(paste("\nConfusion Matrix for", n_genes, "genes:"))
  print(conf_matrix)
  
  # Print results for this gene count
  print(paste("\nResults for", n_genes, "genes:"))
  print(paste("Accuracy:", round(results_df$accuracy[i], 3)))
  print(paste("AUC:", round(results_df$auc[i], 3)))
}

# Print final results table
print("\nFinal Results:")
print(results_df)

# Create plots
# Accuracy and AUC plot
performance_plot <- ggplot(results_df %>% 
                             gather(key = "metric", value = "value", -genes)) +
  geom_line(aes(x = genes, y = value, color = metric)) +
  geom_point(aes(x = genes, y = value, color = metric)) +
  scale_x_log10() +
  labs(title = "KNN Performance vs Number of Genes",
       x = "Number of Genes",
       y = "Value",
       color = "Metric") +
  theme_minimal()

# Save results
ggsave("knn_performance_vs_genes.png", performance_plot)
write.csv(results_df, "knn_binary_classification_results.csv", row.names = FALSE)

# Create summary
summary_text <- paste(
  "Binary Classification (Early Acute vs Others) Results:\n",
  "\nPerformance by gene count:",
  paste("10 genes - Accuracy:", round(results_df$accuracy[1], 3), 
        "AUC:", round(results_df$auc[1], 3)),
  paste("100 genes - Accuracy:", round(results_df$accuracy[2], 3), 
        "AUC:", round(results_df$auc[2], 3)),
  paste("1000 genes - Accuracy:", round(results_df$accuracy[3], 3), 
        "AUC:", round(results_df$auc[3], 3)),
  paste("10000 genes - Accuracy:", round(results_df$accuracy[4], 3), 
        "AUC:", round(results_df$auc[4], 3)),
  sep = "\n"
)

# Print and save summary
cat(summary_text)
writeLines(summary_text, "knn_binary_classification_summary.txt")