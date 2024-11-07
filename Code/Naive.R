# Install and load required packages
install.packages(c("tidymodels", "discrim", "naivebayes", "ggplot2", "reshape2", "data.table"))
library(tidymodels)
library(discrim)  # For naive bayes
library(naivebayes)
library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)

# Load data
expression_data <- fread("HUGO_Final.tsv")
cluster_data <- read.csv("patient_clustering_results.csv")  # Your cluster assignments
metadata <- read.delim("metadata_SRP192714.tsv")  # Your metadata with stages

# Print initial diagnostics
print("Initial data checks:")
print("Expression data dimensions:")
print(dim(expression_data))
print("\nNumber of clusters in Cluster_5000:")
print(table(cluster_data$Cluster_5000))

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

# Set up Naive Bayes model
nb_spec <- 
  naive_Bayes() %>%
  set_mode("classification") %>%
  set_engine("naivebayes")

# Prepare data for cluster prediction
expression_data_t_clusters <- expression_data_t
expression_data_t_clusters$Label <- factor(cluster_data$Cluster_5000[match(rownames(expression_data_t), 
                                                                           cluster_data$Patient_ID)])

# Split data for cluster prediction
set.seed(123)
cluster_split <- initial_split(expression_data_t_clusters, prop = 0.7)
cluster_train <- training(cluster_split)
cluster_test <- testing(cluster_split)

# Train model for cluster prediction
print("\nTraining Naive Bayes model for cluster prediction...")
nb_fit_clusters <- nb_spec %>% 
  fit(Label ~ ., data = cluster_train)

# Make predictions for cluster assignments
all_nb_predictions_cluster <- predict(nb_fit_clusters, expression_data_t)

# Create a new Naive Bayes model for predicting infection stages
# First, create a new data frame with infection stages as labels
expression_data_t_stages <- expression_data_t
expression_data_t_stages$Label <- metadata$refinebio_time[match(rownames(expression_data_t), 
                                                                metadata$refinebio_accession_code)]
expression_data_t_stages$Label <- factor(expression_data_t_stages$Label)

# Train new Naive Bayes model for stages
set.seed(123)
nb_fit_stages <- nb_spec %>% 
  fit(Label ~ ., data = expression_data_t_stages)

# Make predictions for infection stages
all_nb_predictions_stages <- predict(nb_fit_stages, expression_data_t_stages)

# Format predictions for all samples
nb_formatted_all <- data.frame(
  Sample = rownames(expression_data_t),
  Cluster_5000_NB = cluster_assignments,
  Model_5000_NB = as.character(all_nb_predictions_cluster$.pred_class),
  Time_5000_NB = as.character(all_nb_predictions_stages$.pred_class)
)

# Save complete Naive Bayes predictions
write.csv(nb_formatted_all, "nb_predictions_all_samples.csv", row.names = FALSE)

# Make cluster predictions for ALL samples
all_cluster_predictions <- predict(nb_fit_clusters, expression_data_t_clusters)

# Prepare data for stage prediction
expression_data_t_stages <- expression_data_t
expression_data_t_stages$Label <- factor(metadata$refinebio_time[match(rownames(expression_data_t), 
                                                                       metadata$refinebio_accession_code)])

# Train model for stage prediction
print("\nTraining Naive Bayes model for stage prediction...")
nb_fit_stages <- nb_spec %>% 
  fit(Label ~ ., data = expression_data_t_stages)

# Make stage predictions for ALL samples
all_stage_predictions <- predict(nb_fit_stages, expression_data_t_stages)

# Format predictions
nb_formatted_all <- data.frame(
  Sample = rownames(expression_data_t),
  Cluster_5000_NB = cluster_data$Cluster_5000[match(rownames(expression_data_t), 
                                                    cluster_data$Patient_ID)],
  Model_5000_NB = all_cluster_predictions$.pred_class,
  Time_5000_NB = all_stage_predictions$.pred_class
)

# Save formatted predictions
write.csv(nb_formatted_all, "nb_predictions_all_samples.csv", row.names = FALSE)

# Create confusion matrix for cluster predictions
cluster_conf_matrix <- table(Predicted = nb_formatted_all$Model_5000_NB,
                             Actual = nb_formatted_all$Cluster_5000_NB)
print("\nConfusion Matrix for Cluster Predictions:")
print(cluster_conf_matrix)

# Calculate accuracy metrics
cluster_accuracy <- mean(nb_formatted_all$Model_5000_NB == nb_formatted_all$Cluster_5000_NB, 
                         na.rm = TRUE)
print(paste("\nCluster Prediction Accuracy:", round(cluster_accuracy, 3)))

# Print summary statistics
print("\nSummary of Naive Bayes Predictions:")
print("\nCluster Prediction Distribution:")
print(table(nb_formatted_all$Model_5000_NB))
print("\nStage Prediction Distribution:")
print(table(nb_formatted_all$Time_5000_NB))

# Try different numbers of genes
gene_counts <- c(10, 100, 1000, 5000)
accuracies <- data.frame(genes = gene_counts, 
                         cluster_accuracy = NA,
                         stage_accuracy = NA)

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
  
  # Cluster prediction
  current_data_t_clusters <- current_data_t
  current_data_t_clusters$Label <- factor(cluster_data$Cluster_5000[match(rownames(current_data_t), 
                                                                          cluster_data$Patient_ID)])
  
  # Stage prediction
  current_data_t_stages <- current_data_t
  current_data_t_stages$Label <- factor(metadata$refinebio_time[match(rownames(current_data_t), 
                                                                      metadata$refinebio_accession_code)])
  
  # Train and evaluate cluster model
  current_cluster_fit <- nb_spec %>% fit(Label ~ ., data = current_data_t_clusters)
  current_cluster_pred <- predict(current_cluster_fit, current_data_t_clusters)
  cluster_acc <- mean(current_cluster_pred$.pred_class == current_data_t_clusters$Label, 
                      na.rm = TRUE)
  
  # Train and evaluate stage model
  current_stage_fit <- nb_spec %>% fit(Label ~ ., data = current_data_t_stages)
  current_stage_pred <- predict(current_stage_fit, current_data_t_stages)
  stage_acc <- mean(current_stage_pred$.pred_class == current_data_t_stages$Label, 
                    na.rm = TRUE)
  
  # Store results
  accuracies$cluster_accuracy[i] <- cluster_acc
  accuracies$stage_accuracy[i] <- stage_acc
  
  print(paste("Cluster accuracy:", round(cluster_acc, 3)))
  print(paste("Stage accuracy:", round(stage_acc, 3)))
}

# Create accuracy plot
acc_plot <- ggplot(accuracies %>% 
                     gather(key = "prediction_type", value = "accuracy", -genes)) +
  geom_line(aes(x = genes, y = accuracy, color = prediction_type)) +
  geom_point(aes(x = genes, y = accuracy, color = prediction_type)) +
  scale_x_log10() +
  labs(title = "Naive Bayes Model Accuracy vs Number of Genes",
       x = "Number of Genes",
       y = "Accuracy",
       color = "Prediction Type") +
  theme_minimal()

# Save results
write.csv(accuracies, "nb_accuracy_results.csv", row.names = FALSE)
ggsave("nb_accuracy_plot.png", acc_plot)

all_cluster_levels_nb <- unique(c(nb_formatted_all$Model_5000_NB, 
                                  nb_formatted_all$Cluster_5000_NB))
nb_cluster_conf <- confusionMatrix(
  data = factor(nb_formatted_all$Model_5000_NB, levels = all_cluster_levels_nb),
  reference = factor(nb_formatted_all$Cluster_5000_NB, levels = all_cluster_levels_nb)
)
print("\nNaive Bayes Cluster Prediction Confusion Matrix:")
print(nb_cluster_conf)

# For stage predictions
all_stage_levels_nb <- unique(c(nb_formatted_all$Time_5000_NB,
                                metadata$refinebio_time[match(nb_formatted_all$Sample, 
                                                              metadata$refinebio_accession_code)]))
nb_stage_conf <- confusionMatrix(
  data = factor(nb_formatted_all$Time_5000_NB, levels = all_stage_levels_nb),
  reference = factor(metadata$refinebio_time[match(nb_formatted_all$Sample, 
                                                   metadata$refinebio_accession_code)],
                     levels = all_stage_levels_nb)
)
print("\nNaive Bayes Stage Prediction Confusion Matrix:")
print(nb_stage_conf)

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
  
  # Set up and train Naive Bayes model
  nb_spec <- naive_Bayes() %>%
    set_mode("classification") %>%
    set_engine("naivebayes")
  
  # Fit model and predict
  current_fit <- nb_spec %>% fit(Label ~ ., data = current_train)
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
  labs(title = "Naive Bayes Performance vs Number of Genes",
       x = "Number of Genes",
       y = "Value",
       color = "Metric") +
  theme_minimal()

# Save results
ggsave("nb_performance_vs_genes.png", performance_plot)
write.csv(results_df, "nb_binary_classification_results.csv", row.names = FALSE)

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
writeLines(summary_text, "nb_binary_classification_summary.txt")