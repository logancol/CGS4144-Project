library(tidyverse)
library(DESeq2)
library(readr)
library(tibble)
library(pheatmap)
library(data.table)
library(ComplexHeatmap)

#Load and round gene data
df <- readr::read_csv("SRP192714_Hugo.csv")
df_data <- df[3:1023]
df_labels <- df[1:2]
df_rounded <- round(df_data, 0)

#Load and trim metadata
df_meta <- readr::read_tsv("metadata_SRP192714.tsv")
df_meta <- df_meta[c("refinebio_accession_code", "refinebio_time")]
df_meta <- df_meta %>%
  dplyr::mutate(
    refinebio_time = factor(refinebio_time, levels=c("early acute", "late acute", "convalescent"))
  )

#Create DESeq object and run transformation
dds <- DESeqDataSetFromMatrix(
  countData = df_rounded,
  colData = df_meta,
  design = ~refinebio_time
)
dds <- DESeq(dds)
deseq_results <- results(dds)

#Apply log fold change to the results
deseq_results <- lfcShrink(
  dds,
  coef = 3,  # Coefficient for the group comparison of interest
  res = deseq_results,
  type = "normal"  # Use 'normal' shrinkage instead of 'apeglm'
)

#Filter by p-value and logfc
alpha <- 0.05  # Adjusted p-value threshold
log2fc_cutoff <- 1  # Log2 fold change threshold

deseq_results$Hugo <- df_labels[2] 
deseq_results <- as.data.frame(deseq_results) %>%
  rownames_to_column(var="Gene") %>%
  filter(padj < alpha & abs(log2FoldChange) > log2fc_cutoff)

#Generate top 50 genes
top_50_genes <- deseq_results %>%
  arrange(pvalue) %>%
  head(50)

#Set up data for heatmap
heatmap_df <- setDT(data.frame(t(data.frame(assay(dds)))))
filter <- as.numeric(unlist(deseq_results[1]))
heatmap_df <- heatmap_df[,filter,with = FALSE]
heatmap_df <- t(data.frame(heatmap_df))

#Plot heatmap
pheatmap(
  heatmap_df,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df_meta[2],
  main = "Non-Annotated Heatmap",
  scale = "row" # Scale values in the direction of genes (rows)
)