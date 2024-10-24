# Load required libraries
library(tidymodels)
library(broom)
library(ggplot2)
library(dplyr)

# Step 1: Exclude the 'Ensembl' column as it's not numeric and clean up whitespace
top_5000_data_clean <- top_5000_data %>%
  mutate(across(everything(), ~ as.numeric(trimws(.))))  # Convert all columns to numeric

# Step 2: Remove rows with NA, NaN, or Inf values
kmeans_data_clean <- top_5000_data_clean %>%
  filter_all(all_vars(is.finite(.)))

# Step 3: Check if the cleaned data has rows and columns
if (nrow(kmeans_data_clean) == 0 || ncol(kmeans_data_clean) == 0) {
  stop("Error: Cleaned data has no rows or columns.")
} else {
  cat("Cleaned data has", nrow(kmeans_data_clean), "rows and", ncol(kmeans_data_clean), "columns.\n")
}

# Step 4: Perform K-means clustering (k = 3) on cleaned data
set.seed(123)  # For reproducibility
kclust <- kmeans(kmeans_data_clean, centers = 3)

# Step 5: View summary of the clustering results
print(summary(kclust))

# Step 6: Tidy the clustering results
tidy_kmeans <- tidy(kclust)

# Step 7: Add cluster assignments to the original cleaned data
augmented_5000_data <- augment(kclust, kmeans_data_clean)

# Step 8: Run PCA for dimensionality reduction
pca_result <- prcomp(kmeans_data_clean, scale. = TRUE)

# Step 9: Add the first two principal components (PC1 and PC2) to the augmented data
augmented_5000_data$PC1 <- pca_result$x[, 1]
augmented_5000_data$PC2 <- pca_result$x[, 2]

# Step 10: Plot the clusters using PCA results (PC1 and PC2)
ggplot(augmented_5000_data, aes(x = PC1, y = PC2, color = .cluster)) + 
  geom_point() +
  theme_minimal() +
  labs(title = "K-means Clustering of Gene Expression Data (PCA)",
       x = "Principal Component 1",
       y = "Principal Component 2")
