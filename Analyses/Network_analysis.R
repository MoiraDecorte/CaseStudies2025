#loading datasets
efa_data_filtered <- readRDS("data_efa.rds")
efa_data_filtered$Trust_countrymeasure <- as.numeric(efa_data_filtered$Trust_countrymeasure)

# loading packages
library(qgraph)
library(bootnet)
library(EGAnet)
library(dplyr)

## Pairwise Markov Random Field (PMFR) with Trustcountry measure ##

# Make a combined dataset
efa_pmrf_data <- efa_data_filtered %>%
  mutate(Trust_countrymeasure = data_efa$Trust_countrymeasure)

efa_pmrf_data$Trust_countrymeasure <- as.numeric(efa_pmrf_data$Trust_countrymeasure)

pmrf_model <- estimateNetwork(
  efa_pmrf_data,
  default = "EBICglasso",
  corMethod = "cor_auto",  # auto-selects polychoric/Pearson
  tuning = 0.5
)

# Plot the network
plot(pmrf_model,
     layout = "spring",
     title = "PMRF with Trust_countrymeasure",
     labels = colnames(efa_pmrf_data))

# Centrality metrics (strength, betweenness, closeness)
centralityPlot(pmrf_model)

# Bootstrapping for accuracy/stability
pmrftc_boot <- bootnet(
  pmrf_model,
  nBoots = 100,
  type = "case"
)

# Plot edge weight CIs
plot(pmrftc_boot, labels = FALSE, order = "sample")

# Centrality stability (CS-coefficients)
cs_result <- corStability(pmrftc_boot)
print(cs_result)
