library(tidyverse)
library(ggplot2)
library(umap)
library(dplyr)
library(DESeq2)

set.seed(52)

#Load and separate gene data
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

#Norm the data and create the umap table
dds_norm <- vst(dds)
norm_counts <- assay(dds_norm) %>% t()
umap_results <- umap::umap(norm_counts)
umap_plot_df <- data.frame(umap_results$layout) %>%
  tibble::rownames_to_column("refinebio_accession_code") %>%
  dplyr::inner_join(df_meta, by ="refinebio_accession_code")

#Plot the umap data
ggplot(
  umap_plot_df,
  aes(
    x = X1,
    y = X2,
    color = refinebio_time
  )
) +
  geom_point()
