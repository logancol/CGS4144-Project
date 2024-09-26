install.packages("ggplot2")
library(ggplot2)

#Load in csv and log transform data
df <- readr::read_csv("SRP192714_Hugo.csv")
log2_transformed <- mutate_if(df, is.numeric, log2)

#Extract row medians and variances
medians <- c(apply(log2_transformed[3:1023], 1, median))

#plotting density of median log transformed expression values for each gene
ggplot(df, aes(x = medians)) +
  geom_density(fill = "lightblue", color = "darkblue", alpha = 0.6) +
  labs(title = "Density Plot of Median Values", x = "Values", y = "Density") +
  theme_minimal()
