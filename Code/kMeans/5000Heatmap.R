# Assuming you have your data prepared and subsets for top 10, 100, 1000, 5000, and 10000 genes

library(dplyr)
library(tidymodels)
library(ggplot2)

# Perform K-means clustering for each subset
set.seed(123)  # Set seed for reproducibility

# Example clustering for top 10 genes
kclust_10 <- kmeans(top_10_data, centers = 3)  # Adjust 'centers' as needed
cluster_10 <- kclust_10$cluster  # Extract cluster assignments

# Repeat for the other subsets
kclust_100 <- kmeans(top_100_data, centers = 3)
cluster_100 <- kclust_100$cluster

kclust_1000 <- kmeans(top_1000_data, centers = 3)
cluster_1000 <- kclust_1000$cluster

kclust_5000 <- kmeans(top_5000_data, centers = 3)
cluster_5000 <- kclust_5000$cluster

kclust_10000 <- kmeans(top_10000_data, centers = 3)
cluster_10000 <- kclust_10000$cluster

# Combine the clustering results
patient_clustering_results <- data.frame(
  Sample = rownames(top_10_data),  # Assuming the samples (SRR IDs) are the row names
  Cluster_10 = cluster_10,
  Cluster_100 = cluster_100,
  Cluster_1000 = cluster_1000,
  Cluster_5000 = cluster_5000,
  Cluster_10000 = cluster_10000
)

patient_ids_expression <- rownames(patient_clustering_results)
patient_ids_metadata <- metadata$refinebio_accession_code
head(patient_ids_metadata)
head(patient_ids_expression)
metadata_ordered <- metadata[match(patient_ids_expression, metadata$refinebio_accession_code), ]
patient_clustering_results <- cbind(
  metadata_ordered$refinebio_time,
  patient_clustering_results
)

for_heat <- patient_clustering_results[, c("Sample", "Cluster_5000")]
write.table(for_heat, file = "for_heatmap.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# View the organized clustering results
print(patient_clustering_results)
