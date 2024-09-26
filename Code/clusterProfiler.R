if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "DOSE"))

library(clusterProfiler)
library(org.Hs.eg.db)  # Database for human gene annotations
library(DOSE)          
library(tibble)

genes <- read.csv("differentially_expressed_genes.csv")
gene_list_hugo <- as.character(genes$filtered_mapped_hugo_id)

# Convert HUGO gene symbols to Entrez IDs using bitr function
degs_entrez <- bitr(gene_list_hugo, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform GO enrichment analysis
go_enrich <- enrichGO(
  gene = degs_entrez$ENTREZID,     # Use Entrez IDs
  OrgDb = org.Hs.eg.db,            # Database for human genes
  ont = "BP",                      # Ontology: Biological Process (BP)
  pvalueCutoff = 0.05,             # P-value threshold
  pAdjustMethod = "BH",            # Benjamini-Hochberg correction for multiple testing
  qvalueCutoff = 0.05,             # Q-value threshold
  minGSSize = 10,                  # Minimum gene set size
  maxGSSize = 500,                 # Maximum gene set size
  readable = TRUE                  # Converts Entrez IDs back to gene symbols
)

# Explore the enrichment results
head(as.data.frame(go_enrich))

# Visualize the results
dotplot(go_enrich, showCategory = 20)

# Barplot visualization
barplot(go_enrich, showCategory = 20)

# Enrichment map visualization
emapplot(go_enrich)

# Save the GO enrichment results to a CSV file

write.csv(as.data.frame(go_enrich), "GO_enrichment_results.csv")

# Optionally, view a summary of the results in your R environment
summary(go_enrich)
