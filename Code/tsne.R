###TSNE PLOT###

library(readr)
#read meta data
meta <- readr::read_tsv('metadata_SRP192714.tsv')

#get rid of unnecessary columns in meta data
meta <- meta[,-2:-21]
meta <- meta[,-3:-5]

df <- readr::read_csv('SRP192714_Hugo.csv')

#log scale RNA seq data
df_log2 <- df
df_log2[ , -c(1, 2)] <- log2(df_log2[ , -c(1, 2)])
df_new <- df_log2[ , -c(1, 2)]


##install.packages("Rtsne")
#install.packages("ggplot2")

unique_df_new <- unique(df_new)


library(Rtsne)
library(ggplot2)


data <- as.matrix(unique_df_new)  

# Transpose the data
data_t <- t(data)


set.seed(42)  # For testing
tsne_result <- Rtsne(data_t, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 1000)

# use metadata refinebio_time to group results
tsne_df <- data.frame(tsne_result$Y, Time = meta$refinebio_time)
colnames(tsne_df) <- c("Dim1", "Dim2", "Time")

# Plot the t-SNE result with color grouping
ggplot(tsne_df, aes(x = Dim1, y = Dim2, color = Time)) +
  geom_point(size = 2) +
  labs(title = "t-SNE plot colored by Zika infection time") +
  theme_minimal()