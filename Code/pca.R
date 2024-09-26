library(ggplot2)
library(dplyr)
library(DESeq2)

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

#Scale and plot the PCA data
vsd_dds <- varianceStabilizingTransformation(dds, blind = FALSE)
plotPCA(vsd_dds, intgroup = "refinebio_time")