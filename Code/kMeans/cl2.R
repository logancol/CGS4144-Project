# Load required libraries
library(tidymodels)
library(broom)
library(ggplot2)
library(dplyr)

# Assuming 'top_5000_data' or 'top_10000_data' is your data

# Step 1: Convert data to numeric while removing non-numeric entries
# We'll go column by column to handle large datasets and trim whitespaces
top_5000_data_clean <- top_5000_data %>%
  mutate(across(everything(), ~ as.numeric(trimws(.))))

# Step 2: Instead of using filter_all(), we'll create a loop to remove rows with NA, NaN, or Inf efficiently
# The following code removes rows containing NA, NaN, or Inf values

kmeans_data_clean <- top_5000_data_clean[complete.cases(top_5000_data_clean), ]

# Check if there are any rows left after cleaning
if (nrow(kmeans_data_clean) == 0 || ncol(kmeans_data_clean) == 0) {
  stop("Error: Cleaned data has no rows or columns.")
}

# Step 3: Perform K-means clustering (k = 3)
set.seed(123)  # For reproducibility
kclust <- kmeans(kmeans_data_clean, centers = 13)

# Step 4: Augment the data with cluster assignments and PCA for visualization
augmented_5000_data <- augment(kclust, kmeans_data_clean)

# Step 5: Run PCA to reduce the dimensions of the data
pca_result <- prcomp(kmeans_data_clean, scale. = TRUE)

# Step 6: Add the first two principal components (PC1 and PC2) to the augmented data
augmented_5000_data$PC1 <- pca_result$x[, 1]
augmented_5000_data$PC2 <- pca_result$x[, 2]

# Step 7: Plot the clusters using PCA results (PC1 and PC2)
ggplot(augmented_5000_data, aes(x = PC1, y = PC2, color = .cluster)) + 
  geom_point() +
  theme_minimal() +
  labs(title = "K-means Clustering of Gene Expression Data (PCA)",
       x = "Principal Component 1",
       y = "Principal Component 2")
