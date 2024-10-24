

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

# Perform PAM clustering with k clusters (e.g., k = 3)
pam_result <- pam(top_5000_data, k = 3)
pam_result2 <- pam(top_10_data, k = 3)
pam_result3 <- pam(top_100_data, k = 3)
pam_result4 <- pam(top_1000_data, k = 3)
pam_result5 <- pam(genes_10000, k = 3)


pca_result <- prcomp(top_5000_data, scale. = TRUE)
pca_data_5000 <- pca_result$x[, 1:2]
plot(pca_data_5000, col = pam_result$clustering, pch = 19, main = "PAM Clustering 5000 Genes")

pca_result2 <- prcomp(top_10_data, scale. = TRUE)
pca_data_10 <- pca_result2$x[, 1:2]
plot(pca_data_10, col = pam_result2$clustering, pch = 19, main = "PAM Clustering 10 Genes")

pca_result3 <- prcomp(top_100_data, scale. = TRUE)
pca_data_100 <- pca_result3$x[, 1:2]
plot(pca_data_100, col = pam_result3$clustering, pch = 19, main = "PAM Clustering 100 Genes")

pca_result4 <- prcomp(top_1000_data, scale. = TRUE)
pca_data_1000 <- pca_result4$x[, 1:2]
plot(pca_data_1000, col = pam_result4$clustering, pch = 19, main = "PAM Clustering 1000 Genes")

pca_result5 <- prcomp(genes_10000, scale. = TRUE)
pca_data_10000 <- pca_result5$x[, 1:2]
plot(pca_data_10000, col = pam_result5$clustering, pch = 19, main = "PAM Clustering 10000 Genes")

cluster_assignments <- pam_result$clustering
cluster_counts <- table(cluster_assignments)
print(cluster_counts)

# replace the values being passed to c() with the cluster membership numbers for your method
# i.e. if your method finds with 10 genes that cluster 1 has 3 members, Cluster 2 has 4, 
# and cluster 3 has 4, then you should do
# counts_10 <- c(2, 4, 4)
# these are arrays that will represent the results of the clustering method

counts_10 <- c(356,286, 379)
counts_100 <- c(314,319, 388)
counts_1000 <- c(363, 324, 334)
counts_5000 <- c(370, 318, 333)
counts_10000 <- c(309, 353, 359)

# don't change
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

# don't change
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

# creates a data frame to display results, don’t change
results <- data.frame(
  Comparison = c("10 vs 100", "10 vs 1000", "10 vs 5000", "10 vs 10000",
                 "100 vs 1000", "100 vs 5000", "100 vs 10000",
                 "1000 vs 5000", "1000 vs 10000", "5000 vs 10000"),
  P_Value = c(chi_10_100$p.value, chi_10_1000$p.value, chi_10_5000$p.value, chi_10_10000$p.value,
              chi_100_1000$p.value, chi_100_5000$p.value, chi_100_10000$p.value,
              chi_1000_5000$p.value, chi_1000_10000$p.value, chi_5000_10000$p.value)
)
print(results)

#alluvial plot
clusters_10_genes <- pam_result2$clustering
clusters_100_genes <- pam_result3$clustering
clusters_1000_genes <- pam_result4$clustering
clusters_5000_genes <- pam_result$clustering
clusters_10000_genes <- pam_result5$clustering

patient_clustering_results <- data.frame(
  Patient_ID = rownames(t(expression_data)),
  Cluster_10 = clusters_10_genes,
  Cluster_100 = clusters_100_genes,
  Cluster_1000 = clusters_1000_genes,
  Cluster_5000 = clusters_5000_genes,
  Cluster_10000 = clusters_10000_genes
)


alluvial_data <- data.frame(
  Sample = rownames(top_10_data),  # Assuming the samples are the same across datasets
  Cluster_10_genes = as.factor(clusters_10_genes),
  Cluster_100_Genes = as.factor(clusters_100_genes),
  Cluster_1000_Genes = as.factor(clusters_1000_genes),
  Cluster_5000_Genes = as.factor(clusters_5000_genes),
  Cluster_10000_Genes = as.factor(clusters_10000_genes)
)


#install.packages('ggalluvial')
library(ggalluvial)
library(ggplot2)

# Create the alluvial plot
ggplot(alluvial_data, aes(axis1 = Cluster_10_Genes, axis2 = Cluster_100_Genes, 
                          axis3 = Cluster_1000_Genes, axis4 = Cluster_5000_Genes, 
                          axis5 = Cluster_10000_Genes)) +
  geom_alluvium(aes(fill = Cluster_10_Genes)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("10 Genes", "100 Genes", "1000 Genes", "5000 Genes", "10000 Genes"),
                   expand = c(0.1, 0.1)) +
  theme_minimal() +
  labs(fill = "Clusters",
       title = "PAM Clustering: 10, 100, 1000, 5000, and 10000 Genes",
       y = "Number of Samples",
       x = "Gene Set")



#part 4 chi squared test
patient_clustering_results <- data.frame(
  metadata$refinebio_time,
  Sample = rownames(top_10_data),  # Assuming the samples are the same across datasets
  Cluster_10 = as.factor(clusters_10_genes),
  Cluster_100 = as.factor(clusters_100_genes),
  Cluster_1000 = as.factor(clusters_1000_genes),
  Cluster_5000 = as.factor(clusters_5000_genes),
  Cluster_10000 = as.factor(clusters_10000_genes)
)

contingency_10 <- table(patient_clustering_results$Cluster_10, metadata$refinebio_time)
contingency_100 <- table(patient_clustering_results$Cluster_100, metadata$refinebio_time)
contingency_1000 <- table(patient_clustering_results$Cluster_1000, metadata$refinebio_time)
contingency_5000 <- table(patient_clustering_results$Cluster_5000, metadata$refinebio_time)
contingency_10000 <- table(patient_clustering_results$Cluster_10000, metadata$refinebio_time)
chi_test_10 <- chisq.test(contingency_10)
chi_test_100 <- chisq.test(contingency_100)
chi_test_1000 <- chisq.test(contingency_1000)
chi_test_5000 <- chisq.test(contingency_5000)
chi_test_10000 <- chisq.test(contingency_10000)
chi_stage_results <- data.frame(
  Comparison = c("Cluster_10 vs Infection_Stage", "Cluster_100 vs Infection_Stage", 
                 "Cluster_1000 vs Infection_Stage", "Cluster_5000 vs Infection_Stage", 
                 "Cluster_10000 vs Infection_Stage"),
  P_Value = c(chi_test_10$p.value, chi_test_100$p.value, chi_test_1000$p.value, chi_test_5000$p.value, chi_test_10000$p.value)
)
all_results <- rbind(results, chi_stage_results)
adjusted_values <- p.adjust(all_results$P_Value)
all_results$Adjusted_P <- adjusted_values
print(all_results)



