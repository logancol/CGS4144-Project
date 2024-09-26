library(tidyverse)
library(ggplot2)

#Load in csv and log transform data
df <- readr::read_csv("SRP192714_Hugo.csv")
log2_transformed <- mutate_if(df, is.numeric, log2)

#Extract row medians and variances
medians <- c(apply(log2_transformed[3:1023], 1, median))
variances <- c(apply(log2_transformed[3:1023], 1, var))

#Create plot
ggplot(mapping = aes(x = medians, y = variances)) + 
  labs(title="Gene Expression Median vs Gene Expression Variance", x="Gene Expression Median", y="Gene Expression Variance") +
  geom_point()
