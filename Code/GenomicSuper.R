if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GenomicSuperSignature", quietly = TRUE)) {
  BiocManager::install("GenomicSuperSignature", force = TRUE)
}

library(BiocManager)
library(GenomicSuperSignature)


# Loading RAVmodel
# Each model contains an axis of gene expression, so the goal is to see which RAVs contain our
# highly expressed genes.

RAVmodel <- getModel("C2", load = TRUE)

# Prepping df for validation
differentially_expressed_genes <- column_to_rownames(differentially_expressed_genes, var = "Hugo")
expression_matrix <- as.matrix(differentially_expressed_genes)

# Validating differentially expressed genes against RAVmodel data, selecting top 3 RAVs
val_all <- validate(expression_matrix, RAVmodel)
validated_ind <- validatedSignatures(val_all, RAVmodel, num.out = 3, 
                                     swCutoff = 0, indexOnly = TRUE)
enrichment_data <- subsetEnrichedPathways(RAVmodel, validated_ind, include_nes = TRUE)

write.csv(enrichment_data, "GSSEnrichment.csv", row.names = TRUE)

# https://www.bioconductor.org/packages/release/bioc/vignettes/GenomicSuperSignature/inst/doc/Quickstart.html
# followed this example from GSS documentation closely