library(ggplot2)
library(ggalluvial)

# Create a consistent sample subset for each result
common_indices <- 1:10  # Use the same sample size as the smallest clustering result (10 samples)

clusters_10 <- ap_result_10@idx[common_indices]
clusters_100 <- ap_result_100@idx[common_indices]
clusters_1000 <- ap_result_1000@idx[common_indices]
clusters_5000 <- ap_result_5000@idx[common_indices]


# Create a data frame for the alluvial plot
alluvial_data_fixed <- data.frame(
  Sample = common_indices,  # Use consistent sample indices
  Cluster_10 = factor(clusters_10),
  Cluster_100 = factor(clusters_100),
  Cluster_1000 = factor(clusters_1000),
  Cluster_5000 = factor(clusters_5000)
)

# Plot the alluvial diagram
ggplot(alluvial_data_fixed,
       aes(axis1 = Cluster_10, axis2 = Cluster_100, axis3 = Cluster_1000, axis4 = Cluster_5000,
           y = Sample)) +
  geom_alluvium(aes(fill = Cluster_10), width = 0.2) +
  geom_stratum(width = 0.2, fill = "grey", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("10 Genes", "100 Genes", "1000 Genes", "5000 Genes")) +
  theme_minimal() +
  labs(title = "Alluvial Diagram of Cluster Memberships across Gene Subsets",
       x = "Clustering by Gene Subset", y = "Samples") +
  theme(legend.position = "none")

# Save the alluvial plot
ggsave("alluvial_cluster_memberships_fixed.png", width = 12, height = 8)
