contingency_10 <- table(patient_clustering_results$Cluster_10, metadata_ordered$refinebio_time)
contingency_100 <- table(patient_clustering_results$Cluster_100, metadata_ordered$refinebio_time)
contingency_1000 <- table(patient_clustering_results$Cluster_1000, metadata_ordered$refinebio_time)
contingency_5000 <- table(patient_clustering_results$Cluster_5000, metadata_ordered$refinebio_time)
contingency_10000 <- table(patient_clustering_results$Cluster_10000, metadata_ordered$refinebio_time)
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
