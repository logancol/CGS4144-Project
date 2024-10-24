# library init
install.packages("ggalluvial")
install.packages("ggplot2")
library(ggalluvial)
library(ggplot2)
library(reshape2)
library(data.table)
library(tidyclust)
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

# transposing expression data to cluster samples as opposed to individual genes, I didn't just set the df 
# up like this originally because I mistakenly thought that we were clustering genes
expression_data_t <- t(top_10000_data)
colnames(expression_data_t) <- top_10000_genes
column_names <- colnames(expression_data_t)

# subsetting the data to consider different amounts of genes
ten_subset <- expression_data_t[, 1:10] 
hundred_subset <- expression_data_t[, 1:100]  
thousand_subset <- expression_data_t[, 1:1000]  
fivethousand_subset <- expression_data_t[, 1:5000]
tenthousand_subset <- expression_data_t[, 1:10000]  

# clustering implentation https://www.datacamp.com/tutorial/hierarchical-clustering-R
# creating distance matrix for each subset, first step in hierarchical clustering
dm_10 <- dist(ten_subset, method = "euclidean")
dm_100 <- dist(hundred_subset, method = "euclidean")
dm_1000 <- dist(thousand_subset, method = "euclidean")
dm_5000 <- dist(fivethousand_subset, method = "euclidean")
dm_10000 <- dist(tenthousand_subset, method = "euclidean")

# performing clustering for each subset using average linkage method for robustness to noise
hc_10 <- hclust(dm_10, method = "average")
hc_100 <- hclust(dm_100, method = "average")
hc_1000 <- hclust(dm_1000, method = "average")
hc_5000 <- hclust(dm_5000, method = "average")
hc_10000 <- hclust(dm_10000, method = "average")

# cutting the dendrogram (result of hierarchical clustering) to create clusters
clusters_10 <- cutree(hc_10, k = 3)
clusters_100 <- cutree(hc_100, k = 3)
clusters_1000 <- cutree(hc_1000, k = 3)
clusters_5000 <- cutree(hc_5000, k = 3)
clusters_10000 <- cutree(hc_10000, k = 3)

# displaying results
table(clusters_10)
table(clusters_100)
table(clusters_1000)
table(clusters_5000)
table(clusters_10000)

# aggregating cluster membership results for each patient and putting it in a dataframe
patient_clustering_results <- data.frame(
  Sample = rownames(expression_data_t),
  Cluster_10 = clusters_10,
  Cluster_100 = clusters_100,
  Cluster_1000 = clusters_1000,
  Cluster_5000 = clusters_5000,
  Cluster_10000 = clusters_10000
)

# formatting clustering data for each method to be used for alluvial plot
clustering_long <- melt(patient_clustering_results, id.vars = "Sample", variable.name = "Subset_Size", value.name = "Cluster")
clustering_long$Cluster <- as.factor(clustering_long$Cluster) # making cluster number a factor so It can be passed to alluvial plot
# Generate the alluvial plot
ggplot(clustering_long,
       aes(x = Subset_Size, stratum = Cluster, alluvium = Sample, fill = Cluster, label = Cluster)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum() +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  labs(x = "Subset Size", y = "Sample") +
  ggtitle("Alluvial Diagram of Patient Cluster Membership Across Gene Subsets") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(type = "qual", palette = "Set3")


counts_10 <- c(401, 534, 86)
counts_100 <- c(761, 258, 2)
counts_1000 <- c(941, 78, 2)
counts_5000 <- c(933, 86, 2)
counts_10000 <- c(933, 86, 2)
# manually inputting clustering membership results for chi squared test

contingency_10_100 <- rbind(counts_10, counts_100)
contingency_10_1000 <- rbind(counts_10, counts_1000)
contingency_10_5000 <- rbind(counts_10, counts_5000)
contingency_10_10000 <- rbind(counts_10, counts_10000)
contingency_100_1000 <- rbind(counts_100, counts_1000)
contingency_100_5000 <- rbind(counts_100, counts_5000)
contingency_100_10000 <- rbind(counts_100, counts_10000)
contingency_1000_5000 <- rbind(counts_1000, counts_5000)
contingency_1000_10000 <- rbind(counts_1000, counts_10000)
contingency_5000_10000 <- rbind(counts_5000, counts_10000)
# creating contingency matrix for each combination of clustering results for different subset sizes

chi_10_100 <- chisq.test(contingency_10_100)
chi_10_1000 <- chisq.test(contingency_10_1000)
chi_10_5000 <- chisq.test(contingency_10_5000)
chi_10_10000 <- chisq.test(contingency_10_10000)
chi_100_1000 <- chisq.test(contingency_100_1000)
chi_100_5000 <- chisq.test(contingency_100_5000)
chi_100_10000 <- chisq.test(contingency_100_10000)
chi_1000_5000 <- chisq.test(contingency_1000_5000)
chi_1000_10000 <- chisq.test(contingency_1000_10000)
chi_5000_10000 <- chisq.test(contingency_5000_10000)
# performing chi squared test on each contingency matrix

results <- data.frame(
  Comparison = c("10 vs 100", "10 vs 1000", "10 vs 5000", "10 vs 10000",
                 "100 vs 1000", "100 vs 5000", "100 vs 10000",
                 "1000 vs 5000", "1000 vs 10000", "5000 vs 10000"),
  P_Value = c(chi_10_100$p.value, chi_10_1000$p.value, chi_10_5000$p.value, chi_10_10000$p.value,
              chi_100_1000$p.value, chi_100_5000$p.value, chi_100_10000$p.value,
              chi_1000_5000$p.value, chi_1000_10000$p.value, chi_5000_10000$p.value)
)
print(results)
# formatting chi squared results for display

# experimenting with different clustering sizes
clusters_5000_two <- cutree(hc_5000, k = 2)
table(clusters_5000_two)

# statistics

patient_ids_expression <- rownames(expression_data_t)
patient_ids_metadata <- metadata$refinebio_accession_code
head(patient_ids_metadata)
head(patient_ids_expression)
metadata_ordered <- metadata[match(patient_ids_expression, metadata$refinebio_accession_code), ]
# confirming that metadata sample ids and cluster result sample ids are aligned
patient_clustering_results <- cbind(
  metadata_ordered$refinebio_time,
  patient_clustering_results
)
# matching sample clustering assignments with sample refinebio infection times from metadata


contingency_10 <- table(patient_clustering_results$Cluster_10, metadata_ordered$refinebio_time)
contingency_100 <- table(patient_clustering_results$Cluster_100, metadata_ordered$refinebio_time)
contingency_1000 <- table(patient_clustering_results$Cluster_1000, metadata_ordered$refinebio_time)
contingency_5000 <- table(patient_clustering_results$Cluster_5000, metadata_ordered$refinebio_time)
contingency_10000 <- table(patient_clustering_results$Cluster_10000, metadata_ordered$refinebio_time)
# creating contingency tables for cluster membership and infection time 

print(contingency_10)
print(contingency_100)
print(contingency_1000)
print(contingency_5000)
print(contingency_10000)

chi_test_10 <- chisq.test(contingency_10)
chi_test_100 <- chisq.test(contingency_100)
chi_test_1000 <- chisq.test(contingency_1000)
chi_test_5000 <- chisq.test(contingency_5000)
chi_test_10000 <- chisq.test(contingency_10000)
# performing chi squared test on contingency tables
chi_stage_results <- data.frame(
  Comparison = c("Cluster_10 vs Infection_Stage", "Cluster_100 vs Infection_Stage", 
                 "Cluster_1000 vs Infection_Stage", "Cluster_5000 vs Infection_Stage", 
                 "Cluster_10000 vs Infection_Stage"),
  P_Value = c(chi_test_10$p.value, chi_test_100$p.value, chi_test_1000$p.value, chi_test_5000$p.value, chi_test_10000$p.value)
)
all_results <- rbind(results, chi_stage_results)
# formatting chi squared test results comparing infection time to cluser id and combining with the
# results from comparing clustering results over different subset sizes
print(all_results)

adjusted_values <- p.adjust(all_results$P_Value)
all_results$Adjusted_P <- adjusted_values
print(all_results)