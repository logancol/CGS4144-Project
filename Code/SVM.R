install.packages("ggalluvial")
install.packages("ggplot2")
library(ggalluvial)
library(ggplot2)
library(reshape2)
library(data.table)
library(tidyclust)
library(rsample)
library(tidymodels)
library(mclust)
library(dplyr)

# df init

data <- fread("HUGO_Final.tsv")
differentially_expressed <- fread("differentially_expressed_genes.csv")
metadata <- fread("metadata_SRP192714.tsv")

# indexing rows with gene names

rownames(data) <- data$Ensembl
data <- data %>% select(-c(HUGO))
expression_data <- data[, -1, with=FALSE]
rownames(expression_data) <- rownames(data)

# calculating top 10000 genes in variance to be subsetted

variances <- apply(expression_data, 1, var)
sorted_variances <- sort(variances, decreasing = TRUE)
top_10000_indices <- order(variances, decreasing = TRUE)[1:10000]
top_10000_genes <- rownames(expression_data)[top_10000_indices]
top_10000_data <- expression_data[top_10000_indices, ]

# transposing, standardizing, and adding class label (refinebio_time)

expression_data_t <- t(top_10000_data)
colnames(expression_data_t) <- top_10000_genes
patient_ids_metadata <- metadata$refinebio_accession_code
metadata_ordered <- metadata[match(rownames(expression_data_t), metadata$refinebio_accession_code), ]
expression_data_t <- as.data.frame(expression_data_t)
expression_data_t <- as.data.frame(scale(expression_data_t)) # standardize
expression_data_t <- expression_data_t %>%
  mutate(Label = metadata_ordered$refinebio_time)

expression_data_t$Label <- factor(expression_data_t$Label, levels = c("late acute", "convalescent", "early acute"))

ten_subset <- expression_data_t[, c(1:10, 10001)] 
hundred_subset <- expression_data_t[, c(1:100, 10001)]
thousand_subset <- expression_data_t[, c(1:1000, 10001)]
fivethousand_subset <- expression_data_t[, c(1:5000, 10001)]
tenthousand_subset <- expression_data_t[, c(1:10000, 10001)]

# Train SVM model to predict assignment 1 metadata groups

svm_cls_spec <- 
  svm_linear(cost = 1) %>% 
  set_mode("classification") %>% 
  set_engine("kernlab", prob.model = TRUE)
svm_cls_spec


# predicting class labels with 5000 genes 

fivethousand_split <- initial_split(fivethousand_subset, prop=0.7)
fivethousand_training <- training(fivethousand_split)
fivethousand_testing <- testing(fivethousand_split)
set.seed(1)
svm_cls_fit <- svm_cls_spec %>% fit(Label ~ ., data = fivethousand_training)
svm_cls_fit
fivethousand_predicted <- predict(svm_cls_fit, fivethousand_testing) %>% pull(.pred_class)
results <- tibble(
  truth = fivethousand_testing$Label,
  prediction = fivethousand_predicted
)
results %>%
  conf_mat(truth = truth, estimate = prediction)

  
  
# Making class binary so I can do AUC with different gene counts
  
ten_subset$Label <- as.character(ten_subset$Label)
ten_subset$Label <- ifelse(ten_subset$Label %in% c("late acute", "convalescent"), "late acute/convalescent", ten_subset$Label)
ten_subset$Label <- factor(ten_subset$Label)

hundred_subset$Label <- as.character(hundred_subset$Label)
hundred_subset$Label <- ifelse(hundred_subset$Label %in% c("late acute", "convalescent"), "late acute/convalescent", hundred_subset$Label)
hundred_subset$Label <- factor(hundred_subset$Label)

thousand_subset$Label <- as.character(thousand_subset$Label)
thousand_subset$Label <- ifelse(thousand_subset$Label %in% c("late acute", "convalescent"), "late acute/convalescent", thousand_subset$Label)
thousand_subset$Label <- factor(thousand_subset$Label)

tenthousand_subset$Label <- as.character(tenthousand_subset$Label)
tenthousand_subset$Label <- ifelse(tenthousand_subset$Label %in% c("late acute", "convalescent"), "late acute/convalescent", tenthousand_subset$Label)
tenthousand_subset$Label <- factor(tenthousand_subset$Label)

fivethousand_subset_modified <- fivethousand_subset
fivethousand_subset_modified$Label <- as.character(fivethousand_subset_modified$Label)
fivethousand_subset_modified$Label <- ifelse(
  fivethousand_subset_modified$Label %in% c("late acute", "convalescent"),
  "late acute/convalescent",
  fivethousand_subset_modified$Label
)
fivethousand_subset_modified$Label <- factor(fivethousand_subset_modified$Label)

# Predicting class probabilites with 10, 100, 1000, 5000, 10000 genes


# Predicting probabilities with 10 genes
ten_split <- initial_split(ten_subset, prop = 0.7)
ten_training <- training(ten_split)
ten_testing <- testing(ten_split)
svm_cls_fit_ten <- svm_cls_spec_prob %>% fit(Label ~ ., data = ten_training)
ten_predicted_prob <- predict(svm_cls_fit_ten, ten_testing, type = "prob")
results_ten <- ten_testing %>%
  bind_cols(ten_predicted_prob)

# Predicting probabilities with 100 genes
hundred_split <- initial_split(hundred_subset, prop = 0.7)
hundred_training <- training(hundred_split)
hundred_testing <- testing(hundred_split)
svm_cls_fit_hundred <- svm_cls_spec_prob %>% fit(Label ~ ., data = hundred_training)
hundred_predicted_prob <- predict(svm_cls_fit_hundred, hundred_testing, type = "prob")
results_hundred <- hundred_testing %>%
  bind_cols(hundred_predicted_prob)

# Predicting probabilities with 1000 genes
thousand_split <- initial_split(thousand_subset, prop = 0.7)
thousand_training <- training(thousand_split)
thousand_testing <- testing(thousand_split)
svm_cls_fit_thousand <- svm_cls_spec_prob %>% fit(Label ~ ., data = thousand_training)
thousand_predicted_prob <- predict(svm_cls_fit_thousand, thousand_testing, type = "prob")
results_thousand <- thousand_testing %>%
  bind_cols(thousand_predicted_prob)

# Predicting probabilities with 10000 genes
tenthousand_split <- initial_split(tenthousand_subset, prop = 0.7)
tenthousand_training <- training(tenthousand_split)
tenthousand_testing <- testing(tenthousand_split)
svm_cls_fit_tenthousand <- svm_cls_spec_prob %>% fit(Label ~ ., data = tenthousand_training)
tenthousand_predicted_prob <- predict(svm_cls_fit_tenthousand, tenthousand_testing, type = "prob")
results_tenthousand <- tenthousand_testing %>%
bind_cols(tenthousand_predicted_prob)


# Predicting probabilities with 5000 genes (already did labels)
fivethousand_split_modified <- initial_split(fivethousand_subset_modified, prop = 0.7)
fivethousand_training_modified <- training(fivethousand_split_modified)
fivethousand_testing_modified <- testing(fivethousand_split_modified)
svm_cls_fit_fivethousand <- svm_cls_spec_prob %>% fit(Label ~ ., data = fivethousand_training_modified)
fivethousand_predicted_prob <- predict(svm_cls_fit_fivethousand, fivethousand_testing_modified, type = "prob")
results_fivethousand <- fivethousand_testing_modified %>%
bind_cols(fivethousand_predicted_prob)


auc_ten <- results_ten %>%
  roc_auc(truth = Label, ".pred_early acute")
auc_hundred <- results_hundred %>%
  roc_auc(truth = Label, ".pred_early acute")
auc_thousand <- results_thousand %>%
  roc_auc(truth = Label, ".pred_early acute")
auc_fivethousand <- results_fivethousand %>%
  roc_auc(truth = Label, ".pred_early acute")
auc_tenthousand <- results_tenthousand %>%
  roc_auc(truth = Label, ".pred_early acute")

auc_results <- tibble(
  Gene_Count = c(10, 100, 1000, 5000, 10000),
  AUC = c(auc_ten$.estimate, auc_hundred$.estimate, auc_thousand$.estimate, auc_fivethousand$.estimate, auc_tenthousand$.estimate)
)
auc_results



# plotting accuracy comparison between each model
accuracy_fivethousand <- mean(results$truth == results$prediction)
accuracy_ten <- mean(results_ten$truth == results_ten$prediction)
accuracy_hundred <- mean(results_hundred$truth == results_hundred$prediction)
accuracy_thousand <- mean(results_thousand$truth == results_thousand$prediction)
accuracy_tenthousand <- mean(results_tenthousand$truth == results_tenthousand$prediction)
accuracies <- c(accuracy_ten, accuracy_hundred, accuracy_thousand, accuracy_fivethousand, accuracy_tenthousand)
gene_results <- data.frame(
  genes = c(10, 100, 1000, 5000, 10000),
  accuracy = c(accuracy_ten, accuracy_hundred, accuracy_thousand, accuracy_fivethousand, accuracy_tenthousand)
)

acc_plot <- ggplot(gene_results, aes(x = genes, y = accuracy)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  labs(title = "Model Accuracy vs Number of Genes",
       x = "Number of Genes",
       y = "Accuracy") +
  theme_minimal()



cluster_data <- fread("for_heatmap.tsv")

# predicting clusters from assignment 3 with 5000 genes

fivethousand_subset$Label <- factor(cluster_data$Cluster_5000, labels = c(1, 2, 3))
fivethousand_split <- initial_split(fivethousand_subset, prop=0.7)
fivethousand_training <- training(fivethousand_split)
fivethousand_testing <- testing(fivethousand_split)
set.seed(1)
svm_cls_fit <- svm_cls_spec %>% fit(Label ~ ., data = fivethousand_training)
svm_cls_fit
fivethousand_predicted <- predict(svm_cls_fit, fivethousand_testing) %>% pull(.pred_class)


# writing class label and cluster label results to file to compare with ROC curve

sample_names <- rownames(fivethousand_subset)
predictive_results <- data.frame(
  Sample = sample_names,
  Cluster_5000_SVM = fivethousand_predicted,
  Model_5000_SVM = fivethousand_predicted_total_model
)
write.table(predictive_results, "predictive_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

