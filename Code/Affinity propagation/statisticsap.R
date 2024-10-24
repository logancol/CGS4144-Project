samples <- metadata$refinebio_accession_code

clusters_10 <- ap_result_10@idx
clusters_100 <- ap_result_100@idx
clusters_1000 <- ap_result_1000@idx
clusters_5000 <- ap_result_5000@idx

clusters_10 <- clusters_10[1:length(samples)]
clusters_100 <- clusters_100[1:length(samples)]
clusters_1000 <- clusters_1000[1:length(samples)]
clusters_5000 <- clusters_5000[1:length(samples)]

cluster_df <- data.frame(
  Sample = samples, 
  Group = metadata$refinebio_time,
  Cluster_10 = clusters_10,
  Cluster_100 = clusters_100,
  Cluster_1000 = clusters_1000,
  Cluster_5000 = clusters_5000
)

head(cluster_df)
View(cluster_df)
write.csv(cluster_df, "cluster_memberships.csv", row.names = FALSE)

patient_clustering_results <- cluster_df

contingency_10 <- table(patient_clustering_results$Cluster_10, metadata_ordered$refinebio_time)
contingency_100 <- table(patient_clustering_results$Cluster_100, metadata_ordered$refinebio_time)
contingency_1000 <- table(patient_clustering_results$Cluster_1000, metadata_ordered$refinebio_time)
contingency_5000 <- table(patient_clustering_results$Cluster_5000, metadata_ordered$refinebio_time)

chi_test_10 <- chisq.test(contingency_10)
chi_test_100 <- chisq.test(contingency_100)
chi_test_1000 <- chisq.test(contingency_1000)
chi_test_5000 <- chisq.test(contingency_5000)

chi_stage_results <- data.frame(
  Comparison = c("Cluster_10 vs Infection_Stage", "Cluster_100 vs Infection_Stage", 
                 "Cluster_1000 vs Infection_Stage", "Cluster_5000 vs Infection_Stage"),
  P_Value = c(chi_test_10$p.value, chi_test_100$p.value, chi_test_1000$p.value, chi_test_5000$p.value)
)

all_results <- rbind(results, chi_stage_results)
adjusted_values <- p.adjust(all_results$P_Value)
all_results$Adjusted_P <- adjusted_values

print(all_results)
