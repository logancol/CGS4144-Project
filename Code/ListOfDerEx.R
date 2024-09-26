# Load necessary libraries for data manipulation and differential expression analysis
library(tidyverse)  # For data wrangling
library(DESeq2)     # For differential expression analysis
library(readr)      # For reading CSV and TSV files
library(tibble)     # For tibble manipulation

# Load the gene expression data from CSV file
df <- readr::read_csv("SRP192714_Hugo.csv")

# Remove unnecessary column 'all_hugo_ids' from the dataset
df <- subset(df, select = -all_hugo_ids)  

# Round all numeric columns in the dataset to whole numbers (assuming counts)
df_rounded <- df %>% mutate_if(is.numeric, round)  

# Extract the expression data for differential analysis (assuming columns 3-1023 contain expression data)
df_rounded_filtered <- df_rounded[3:1023]  

# Extract the first two columns that likely contain gene names and other metadata
df_gene_names <- df_rounded[1:2]

# Load the metadata that contains sample information (like condition/time)
df_meta <- readr::read_tsv("metadata_SRP192714.tsv")

# Select only relevant columns from metadata (sample ID and time point) and order the factor levels for the time variable
df_meta <- df_meta %>%
  select(refinebio_accession_code, refinebio_time) %>%  
  mutate(refinebio_time = factor(refinebio_time, levels = c("early acute", "late acute", "convalescent")))

# Rename the columns of the filtered expression data to match the sample IDs from metadata
colnames(df_rounded_filtered) <- df_meta$refinebio_accession_code

# Create a DESeq2 dataset object using the expression data and metadata
ddset <- DESeqDataSetFromMatrix(
  countData = as.matrix(df_rounded_filtered),  # Expression data matrix (count data)
  colData = df_meta,                          # Metadata for the samples
  design = ~ refinebio_time                   # Design formula to model the effect of 'refinebio_time' (time points)
)

# Perform differential expression analysis using DESeq2
deseq_object <- DESeq(ddset)

# Extract the results of the differential expression analysis
deseq_results <- results(deseq_object)

# Shrink the log2 fold changes to avoid inflated values for genes with low counts
deseq_results_other <- lfcShrink(
  deseq_object,
  coef = 3, 
  res = deseq_results,
  type = "normal"  # Use normal shrinkage method
)

# Add Hugo gene names back to the results for easier interpretation
deseq_results$Hugo <- df_gene_names[2]

# Set the significance thresholds for adjusted p-value and log2 fold change
alpha <- 0.05               
log2fc_cutoff <- 1          

# Filter the DESeq2 results to include only significant genes that pass the p-value and fold change thresholds
deseq_results <- as.data.frame(deseq_results) %>%
  rownames_to_column(var = "Gene") %>%  # Convert row names (gene IDs) to a column for easier manipulation
  filter(padj < alpha & abs(log2FoldChange) > log2fc_cutoff)  # Apply significance filters

# Show the first few rows of the filtered, significant results
head(deseq_results)

# Write the differentially expressed genes to a CSV file
write.csv(deseq_results, "differentially_expressed_genes.csv", row.names = TRUE)
