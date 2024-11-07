library(reshape2)
library(ggalluvial)
library(ggplot2)
library(reshape2)
library(data.table)
library(tidyclust)
library(rsample)
library(tidymodels)
library(mclust)
library(dplyr)
if (!require(pROC)) install.packages("pROC")
library(pROC)

# reading and formatting data from each model results

SVM_results <- fread("predictive_results.tsv")
KNN_results <- fread("knn_predictions_all_samples.csv")
KNN_results$Model_5000_KNN <- NULL
colnames(KNN_results)[3] <- "Model_5000_KNN"
Naive_results <- fread("nb_predictions_all_samples.csv")
Naive_results$Cluster_5000_NB <- NULL
colnames(Naive_results)[2] <- "Cluster_5000_NB"
colnames(Naive_results)[3] <- "Model_5000_NB"
LOG_results <- fread("data_carmen.csv")
LOG_results$V1 <- NULL
colnames(LOG_results)[1] <- "Sample"
colnames(LOG_results)[2] <- "Cluster_5000_LOG"
colnames(LOG_results)[3] <- "Model_5000_LOG"
RF_results <- fread("rf_predictions_all_samples.csv")
RF_results$Model_5000_RF <- NULL
colnames(RF_results)[3] <- "Model_5000_RF"


# formatting df containing results from each model

main_results <- merge(SVM_results, KNN_results, by = "Sample")
main_results <- merge(main_results, Naive_results, by = "Sample")
main_results <- merge(LOG_results, main_results, by = "Sample")
main_results <- merge(RF_results, main_results, by = "Sample")
metadata <- fread("metadata_SRP192714.tsv")
colnames(metadata)[1] <- "Sample"
metadata_sub <- metadata[, c("Sample", "refinebio_time")]
main_results <- merge(main_results, metadata_sub, by = "Sample")

# melting class data

class_long <- melt(main_results, id.vars = "Sample", measure.vars = grep("Model_", names(main_results), value = TRUE), variable.name = "Model", value.name = "Class_Label")
class_label_counts <- as.data.frame.matrix(table(class_long$Sample, class_long$Class_Label))
class_label_counts <- cbind(Sample = rownames(class_label_counts), class_label_counts)
rownames(class_label_counts) <- NULL

# melting clustering data

cluster_long <- melt(main_results, id.vars = "Sample", measure.vars = grep("Cluster_", names(main_results), value = TRUE), variable.name = "Model", value.name = "Cluster")
cluster_counts <- as.data.frame.matrix(table(cluster_long$Sample, cluster_long$Cluster))
cluster_counts <- cbind(Sample = rownames(cluster_counts), cluster_counts)
rownames(cluster_counts) <- NULL

# finding stability as amount of models that predicted correct label divided by total models

class_stability <- class_label_counts
class_stability$Class_Stability <- apply(class_label_counts[, -1], 1, function(x) max(x) / sum(x))

cluster_stability <- cluster_counts
cluster_stability$Cluster_Stability <- apply(cluster_counts[, -1], 1, function(x) max(x) / sum(x))

stability_df <- class_stability %>%
  select(Sample, Class_Stability) %>%
  left_join(select(cluster_stability, Sample, Cluster_Stability), by = "Sample")

# correlation test with multiple test correction

p_values <- c()
p_values <- c(p_values, cor.test(stability_df$Class_Stability, stability_df$Cluster_Stability)$p.value)
padj_values <- p.adjust(p_values, method = "BH")
print(padj_values)

# padj value lower than 0.05, we can infer some correlation


# Combining two of the labels (late acute and convalecent) into one for AUC calculation

main_results_AUC <- main_results %>%
  mutate(refinebio_time = ifelse(refinebio_time %in% c("late acute", "convalescent"), 
                                 "late acute/convalescent", 
                                 refinebio_time))
main_results_AUC <- main_results_AUC %>%
  mutate(Model_5000_KNN = ifelse(Model_5000_KNN %in% c("late acute", "convalescent"), 
                                 "late acute/convalescent", 
                                 Model_5000_KNN))

main_results_AUC <- main_results_AUC %>%
  mutate(Model_5000_NB = ifelse(Model_5000_NB %in% c("late acute", "convalescent"), 
                                 "late acute/convalescent", 
                                 Model_5000_NB))

main_results_AUC <- main_results_AUC %>%
  mutate(Model_5000_LOG = ifelse(Model_5000_LOG %in% c("late acute", "convalescent"), 
                                 "late acute/convalescent", 
                                 Model_5000_LOG))

main_results_AUC <- main_results_AUC %>%
  mutate(Model_5000_SVM = ifelse(Model_5000_SVM %in% c("late acute", "convalescent"), 
                                 "late acute/convalescent", 
                                 Model_5000_SVM))

main_results_AUC <- main_results_AUC %>%
  mutate(Model_5000_RF = ifelse(Model_5000_RF %in% c("late acute", "convalescent"), 
                                 "late acute/convalescent", 
                                 Model_5000_RF))


