#loading packages
library(cluster)
library(factoextra)
library(dplyr)

#loading datasets
factor_scores <- readRDS("factor_scores.rds")

# Select columns: MR1–MR9, Trust, Continent, etc.
factor_scores_clean <- factor_scores %>%
  select(starts_with("MR"), Trust_countrymeasure, Continent, SubjectID, country) %>%
  na.omit()

colnames(factor_scores_clean)

# Scale factor scores
scaled_data <- scale(factor_scores_clean[, grep("^MR", colnames(factor_scores_clean))])

# determine optimal number of clusters (elbow or silhouette)
fviz_nbclust(scaled_data, kmeans, method = "silhouette")  # or method = "wss"

# Perform k-means clustering with your chosen number of clusters (e.g., 2)
set.seed(123)
kmeans_result <- kmeans(scaled_data, centers = 2, nstart = 25, iter.max = 100)
# Add cluster assignments to your data
factor_scores_clean$cluster <- factor(kmeans_result$cluster)


# Visualize clusters
fviz_cluster(kmeans_result, data = scaled_data,
             palette = "jco", ggtheme = theme_minimal(),
             main = "K-means Clustering on Psychological Profiles")

# Summarize cluster profiles
cluster_profiles <- factor_scores_clean %>%
  group_by(cluster) %>%
  summarize(across(starts_with("MR"), list(mean = mean, sd = sd)))

View(cluster_profiles)

# Summarize clusters per continent
factor_scores_clean %>%
  count(cluster, Continent) %>%
  tidyr::pivot_wider(names_from = Continent, values_from = n, values_fill = 0)

factor_scores_clean %>%
  group_by(cluster, Continent) %>%
  summarise(n = n()) %>%
  group_by(cluster) %>%
  mutate(percent = round(n / sum(n) * 100, 1)) %>%
  arrange(cluster)

# Calculate mean Trust_countrymeasure per cluster
mean_trust_by_cluster <- factor_scores_clean %>%
  group_by(cluster) %>%
  summarize(
    mean_trust = mean(Trust_countrymeasure, na.rm = TRUE),
    sd_trust = sd(Trust_countrymeasure, na.rm = TRUE),
    n = n()
  )

print(mean_trust_by_cluster)

