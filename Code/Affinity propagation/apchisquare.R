# Extracting cluster assignments and counts for each result

# For 10 genes
ap_clusters_10 <- ap_result_10@idx  # Get cluster assignments for 10 genes
cluster_counts_10 <- table(ap_clusters_10)  # Count genes per cluster
counts_10 <- as.vector(cluster_counts_10)

# For 100 genes
ap_clusters_100 <- ap_result_100@idx  # Get cluster assignments for 100 genes
cluster_counts_100 <- table(ap_clusters_100)  # Count genes per cluster
counts_100 <- as.vector(cluster_counts_100)

# For 1000 genes
ap_clusters_1000 <- ap_result_1000@idx  # Get cluster assignments for 1000 genes
cluster_counts_1000 <- table(ap_clusters_1000)  # Count genes per cluster
counts_1000 <- as.vector(cluster_counts_1000)

# For 5000 genes
ap_clusters_5000 <- ap_result_5000@idx  # Get cluster assignments for 5000 genes
cluster_counts_5000 <- table(ap_clusters_5000)  # Count genes per cluster
counts_5000 <- as.vector(cluster_counts_5000)

# Pad counts if the number of clusters differs
max_clusters <- max(length(counts_10), length(counts_100), length(counts_1000), length(counts_5000))
counts_10 <- c(counts_10, rep(0, max_clusters - length(counts_10)))
counts_100 <- c(counts_100, rep(0, max_clusters - length(counts_100)))
counts_1000 <- c(counts_1000, rep(0, max_clusters - length(counts_1000)))
counts_5000 <- c(counts_5000, rep(0, max_clusters - length(counts_5000)))

# Running chi-squared tests on each pair of clustering results

contingency_10_100 <- rbind(counts_10, counts_100)
contingency_10_1000 <- rbind(counts_10, counts_1000)
contingency_10_5000 <- rbind(counts_10, counts_5000)
contingency_100_1000 <- rbind(counts_100, counts_1000)
contingency_100_5000 <- rbind(counts_100, counts_5000)
contingency_1000_5000 <- rbind(counts_1000, counts_5000)

chi_10_100 <- chisq.test(contingency_10_100)
chi_10_1000 <- chisq.test(contingency_10_1000)
chi_10_5000 <- chisq.test(contingency_10_5000)
chi_100_1000 <- chisq.test(contingency_100_1000)
chi_100_5000 <- chisq.test(contingency_100_5000)
chi_1000_5000 <- chisq.test(contingency_1000_5000)

# Create a data frame to display results
results <- data.frame(
  Comparison = c("10 vs 100", "10 vs 1000", "10 vs 5000", 
                 "100 vs 1000", "100 vs 5000", "1000 vs 5000"),
  P_Value = c(chi_10_100$p.value, chi_10_1000$p.value, chi_10_5000$p.value,
              chi_100_1000$p.value, chi_100_5000$p.value, chi_1000_5000$p.value)
)
print(results)
