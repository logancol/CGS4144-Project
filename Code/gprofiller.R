# Load necessary library
library(gprofiler2)

# Load the list of differentially expressed genes from the CSV
genes <- read.csv("differentially_expressed_genes.csv")

# Convert gene IDs to character vector
gene_list <- as.character(genes$filtered_mapped_hugo_id)

# Perform functional enrichment analysis using gost
gost_results <- gost(
  query = gene_list,              # Your list of genes
  organism = "hsapiens",          # Set to "hsapiens" for human genes
  sources = c("GO:BP"),   # Limit to GO:BP and Reactome
  ordered_query = FALSE,          # Set to TRUE if the list is ordered by importance
  multi_query = FALSE,            # Single query for one list of genes
  significant = TRUE,            # Show all results, including non-significant
  exclude_iea = FALSE,            # Include in-silico annotations
  user_threshold = 0.05,          # P-value threshold
  correction_method = "g_SCS"     # Multiple testing correction (default g:SCS)
)

# Extract the result from gost_results
enrichment_table <- gost_results$result

# Check for list columns and flatten them (especially the "intersection" column)
enrichment_table[] <- lapply(enrichment_table, function(x) {
  if (is.list(x)) {
    sapply(x, function(y) paste(y, collapse = ", "))
  } else {
    x
  }
})

# Create a Manhattan plot of the enrichment results
gostplot(gost_results, capped = TRUE, interactive = TRUE)


# Save the enrichment table as a CSV file
write.csv(enrichment_table, "gprofiler_enrichment_results.csv", row.names = FALSE)
