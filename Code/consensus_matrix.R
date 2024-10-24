library(tidyclust)
library(tidyverse)
library(DESeq2)
library(BiocManager)
library(ConsensusClusterPlus)
library(ggalluvial)
library(reshape)
library(pheatmap)
library(scales)
library(RColorBrewer)
library(ComplexHeatmap)

#Load and round gene data
df <- readr::read_csv("SRP192714_Hugo.csv")
df_data <- df[3:1023]
df_labels <- df[1:2]

rownames(df) <- df$Ensembl
variances <- apply(df_data, 1, var)

indexes <- order(variances, decreasing = TRUE)[1:10000]
top_10000_genes <- rownames(df)[indexes]
top_10000 <- df[rownames(df) %in% top_10000_genes,]

indexes <- order(variances, decreasing = TRUE)[1:5000]
top_5000_genes <- rownames(df)[indexes]
top_5000 <- df[rownames(df) %in% top_5000_genes,]

indexes <- order(variances, decreasing = TRUE)[1:1000]
top_1000_genes <- rownames(df)[indexes]
top_1000 <- df[rownames(df) %in% top_1000_genes,]

indexes <- order(variances, decreasing = TRUE)[1:100]
top_100_genes <- rownames(df)[indexes]
top_100 <- df[rownames(df) %in% top_100_genes,]

indexes <- order(variances, decreasing = TRUE)[1:10]
top_10_genes <- rownames(df)[indexes]
top_10 <- df[rownames(df) %in% top_10_genes,]



  cluster10000 <- ConsensusClusterPlus(
    data.matrix(top_10000[3:1023]), maxK=6, reps=100, pItem=0.8, pFeature=1, plot=NULL)
  
  cluster5000 <- ConsensusClusterPlus(
    data.matrix(top_5000[3:1023]), maxK=6, reps=100, pItem=0.8, pFeature=1, plot=NULL)
  
  cluster1000 <- ConsensusClusterPlus(
    data.matrix(top_1000[3:1023]), maxK=6, reps=100, pItem=0.8, pFeature=1, plot=NULL)
  
  cluster100 <- ConsensusClusterPlus(
    data.matrix(top_100[3:1023]), maxK=6, reps=100, pItem=0.8, pFeature=1, plot=NULL)
  
  cluster10 <- ConsensusClusterPlus(
    data.matrix(top_10[3:1023]), maxK=6, reps=100, pItem=0.8, pFeature=1, plot=NULL)

counts_10000 <- table(cluster10000[[3]]$consensusClass)
counts_5000 <- table(cluster5000[[3]]$consensusClass)
counts_1000 <- table(cluster1000[[3]]$consensusClass)
counts_100 <- table(cluster100[[3]]$consensusClass)
counts_10 <- table(cluster10[[3]]$consensusClass)

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

results <- data.frame(
  Comparison = c("10 vs 100", "10 vs 1000", "10 vs 5000", "10 vs 10000",
                 "100 vs 1000", "100 vs 5000", "100 vs 10000",
                 "1000 vs 5000", "1000 vs 10000", "5000 vs 10000"),
  P_Value = c(chi_10_100$p.value, chi_10_1000$p.value, chi_10_5000$p.value, chi_10_10000$p.value,
              chi_100_1000$p.value, chi_100_5000$p.value, chi_100_10000$p.value,
              chi_1000_5000$p.value, chi_1000_10000$p.value, chi_5000_10000$p.value)
)

results$P_Value <- p.adjust(results$P_Value)


clustering_results <- data.frame(
  Patient_ID = colnames(df_data),
  Cluster10 <- cluster10[[3]]$consensusClass,
  Cluster100 <- cluster100[[3]]$consensusClass,
  Cluster1000 <- cluster1000[[3]]$consensusClass,
  Cluster5000 <- cluster5000[[3]]$consensusClass,
  Cluster10000 <- cluster10000[[3]]$consensusClass
)

colnames(clustering_results) <- c("Patient_ID", "Cluster10", "Cluster100", 
                                  "Cluster1000", "Cluster5000", "Cluster10000")

melted_clustering <- melt(clustering_results, id.vars="Patient_ID")
colnames(melted_clustering) <- c("Patient_ID", "Size", "Cluster")

melted_clustering$Cluster <- as.factor(melted_clustering$Cluster)

ggplot(melted_clustering,
       aes(x=Size, stratum=Cluster, alluvium=Patient_ID, fill=Cluster, label=Cluster)) +
  geom_flow(stat="alluvium", lode.guidance="frontback", color="darkgray")+
  geom_stratum()+
  theme_minimal()+
  ggtitle("Alluvial Diagram of Cluster Membership Across Gene Subsets with Consensus Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_brewer(type = "qual", palette = "Set3")

df_meta <- readr::read_tsv("metadata_SRP192714.tsv")
df_meta <- df_meta[c("refinebio_accession_code", "refinebio_time")]
df_meta <- df_meta %>%
  dplyr::mutate(
    refinebio_time = factor(refinebio_time, levels=c("early acute", "late acute", "convalescent"))
  )

clustering_results$refinebio_time <- df_meta$refinebio_time

contingency_10000_time <- table(clustering_results$Cluster10000, clustering_results$refinebio_time)
contingency_5000_time <- table(clustering_results$Cluster5000, clustering_results$refinebio_time)
contingency_1000_time <- table(clustering_results$Cluster1000, clustering_results$refinebio_time)
contingency_100_time <- table(clustering_results$Cluster100, clustering_results$refinebio_time)
contingency_10_time <- table(clustering_results$Cluster10, clustering_results$refinebio_time)

chi_10000_time <- chisq.test(contingency_10000_time)
chi_5000_time <- chisq.test(contingency_5000_time)
chi_1000_time <- chisq.test(contingency_1000_time)
chi_100_time <- chisq.test(contingency_100_time)
chi_10_time <- chisq.test(contingency_10_time)

results <- data.frame(
  Comparison = c("10 vs Time", "100 vs Time", "1000 vs Time", "5000 vs Time", "1000 vs Time"),
  P_Value = c(chi_10_time$p.value, chi_100_time$p.value, chi_1000_time$p.value, 
              chi_5000_time$p.value, chi_10000_time$p.value)
)

results$P_Value <- p.adjust(results$P_Value)

k_means_cluster <- readr::read_tsv("k_means_cluster.tsv")
h_cluster <- readr::read_tsv("h_cluster.tsv")
pam_cluster <- readr::read_csv("pam_clustering.csv")
ap_cluster <- readr::read_tsv("apcluster.tsv")

cluster5000[[3]]$consensusClass <- as.factor(cluster5000[[3]]$consensusClass)
k_means_cluster$Cluster_5000 <- as.factor(k_means_cluster$Cluster_5000)
h_cluster$Cluster_5000 <- as.factor(h_cluster$Cluster_5000)
pam_cluster$Cluster_5000 <- as.factor(pam_cluster$Cluster_5000)
ap_cluster$Cluster_5000 <- as.factor(ap_cluster$Cluster_5000)

annotations <- data.frame(
  df_meta$refinebio_time,
  cluster5000[[3]]$consensusClass,
  k_means_cluster$Cluster_5000,
  h_cluster$Cluster_5000,
  pam_cluster$Cluster_5000,
  ap_cluster$Cluster_5000
)

colnames(annotations) <- c("Infection Stage", "Consensus",
                           "Kmeans", "HCluster", "PAM", "AffinityPropagation")


ComplexHeatmap::Heatmap(
  as.matrix(top_5000[3:1023]),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_title = "Samples",
  row_title = "Genes",
  top_annotation  = HeatmapAnnotation(df = annotations, annotation_name_side ="left"),
  use_raster = TRUE,
  raster_quality = 5,
  heatmap_width = unit(8.5, "cm")
)
