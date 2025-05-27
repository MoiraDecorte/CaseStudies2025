#Packages laden
library(readr)
library(dplyr)
library(psych)
library(GGally)
library(ggplot2)
library(gridExtra)

#Data inlezen
data_raw <- read_csv( "/Users/ann-sophie/Desktop/Case studies 2025/COVIDiSTRESS global survey May 30 2020 (___final cleaned file___).csv",
  locale = locale(encoding = "UTF-8"), 
  na = c("NA", "")
)

#Kolomnamen netjes maken
spaces_to_underscores <- function(data) {
  names(data) <- gsub("[ -]", "_", names(data))
  return(data)
}

data_clean <- spaces_to_underscores(data_raw)

#Alleen Likert-schalen selecteren (1–5 en 1–6)
likert_items <- grep("^(Scale_PSS10_UCLA_|Corona_concerns_|Compliance_|BFF_15_|Expl_Distress_|SPS_|Expl_Coping_|Expl_media_)", 
                     names(data_clean), value = TRUE)

oecd_items <- grep("^(OECD_people_|OECD_insititutions_)", names(data_clean), value = TRUE)

# Voeg Trust_countrymeasure expliciet toe
extra_items <- c("Trust_countrymeasure")

all_items <- c(likert_items, oecd_items, extra_items)

efa_data <- data_clean[, all_items]

distress_vars <- grep("^Expl_Distress_", names(efa_data), value = TRUE)

efa_data[distress_vars] <- lapply(efa_data[distress_vars], function(x) {
  x[x == 7] <- NA
  return(x)
})

efa_data[efa_data == 99] <- NA

# Alles numeriek maken
efa_data <- efa_data[ , names(efa_data) != "Expl_Distress_txt" ]
efa_data[] <- lapply(efa_data, function(x) as.numeric(as.character(x)))

# STANDAARDISEREN: zet alles op z-score
efa_data_scaled <- as.data.frame(scale(efa_data))

non_numeric_check <- sapply(efa_data, function(x) {
  any(!grepl("^[0-9\\.]+$", na.omit(as.character(x))))
})
names(efa_data)[non_numeric_check]
na_inducing_values <- lapply(efa_data, function(x) {
  unique(x[is.na(suppressWarnings(as.numeric(as.character(x)))) & !is.na(x)])
})
na_inducing_values <- Filter(length, na_inducing_values)
print(na_inducing_values)


#Verwijder variabelen zonder voldoende correlatie
remove_uncorrelated_vars <- function(df, threshold = 0.3) {
  cor_matrix <- cor(df, use = "pairwise.complete.obs")
  has_correlations <- sapply(1:ncol(cor_matrix), function(i) {
    any(abs(cor_matrix[i, -i]) >= threshold)
  })
  df_filtered <- df[, has_correlations]
  return(df_filtered)
}
efa_data_filtered <- remove_uncorrelated_vars(efa_data_scaled, threshold = 0.3)


# Flexibele correlatiematrixvisualisatie (detecteert of z-scores of ruwe Likert-data)
ggplotMatrix <- function(efa_data, bin_width = 0.5) {
  # Slimme detectie van schaalbereik
  data_range <- range(efa_data, na.rm = TRUE)
  if (data_range[2] > 6) {
    axis_range <- c(0, 10)  # voor brede schalen
  } else if (data_range[1] < 0) {
    axis_range <- c(-3, 3)  # voor gestandaardiseerde data
  } else {
    axis_range <- c(1, 6)   # klassieke Likert
  }
  my_bins <- seq(axis_range[1], axis_range[2], bin_width)
  axis_breaks <- floor(axis_range[1]):ceiling(axis_range[2])
  
  my_heatmap <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      stat_bin2d(breaks = list(x = my_bins, y = my_bins)) +
      coord_cartesian(xlim = axis_range, ylim = axis_range) +
      scale_x_continuous(breaks = axis_breaks) +
      scale_y_continuous(breaks = axis_breaks) +
      scale_fill_gradient(low = "white", high = "red") +
      theme(legend.position = "none")
  }
  
  my_diag <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_histogram(breaks = my_bins) +
      coord_cartesian(xlim = axis_range) +
      scale_x_continuous(breaks = axis_breaks)
  }
  
  my_upper <- function(data, mapping, ...) {
    x <- eval_data_col(data, mapping$x)
    y <- eval_data_col(data, mapping$y)
    corr <- cor(x, y, use = "pairwise.complete.obs")
    
    ggplot(data = data, mapping = mapping) +
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
                fill = colorRampPalette(c("blue", "white", "red"))(100)[
                  findInterval(corr, seq(-1, 1, length.out = 100))
                ]) +
      annotate("text", x = mean(range(x, na.rm = TRUE)),
               y = mean(range(y, na.rm = TRUE)),
               label = sprintf("%.2f", corr),
               size = 3.5) +
      theme_void()
  }
  
  ggpairs(efa_data,
          lower = list(continuous = my_heatmap),
          diag = list(continuous = my_diag),
          upper = list(continuous = my_upper)) +
    theme_bw()
}

# Correlatiematrix voor efa_data_filtered
cor_matrix <- cor(efa_data_filtered, use = "pairwise.complete.obs")
heatmap(cor_matrix)

# Mooie correlatieplot
library(corrplot)
corrplot(cor_matrix, method = "color", tl.cex = 0.5)

# Overzicht van verwijderde variabelen (te lage correlaties)
removed_vars <- setdiff(names(efa_data_scaled), names(efa_data_filtered))

# Toon waarom ze zijn verwijderd
removed_vars_numeric <- intersect(removed_vars, colnames(efa_data_scaled))
cor_matrix_all <- cor(efa_data_scaled, use = "pairwise.complete.obs")

cor_details <- data.frame(
  variable = removed_vars_numeric,
  max_correlation = sapply(removed_vars_numeric, function(var) {
    max(abs(cor_matrix_all[var, setdiff(colnames(cor_matrix_all), var)]), na.rm = TRUE)
  })
)

# Sorteer en toon de zwakst gecorreleerde eerst
cor_details <- cor_details[order(cor_details$max_correlation), ]
print(cor_details)

#  Matrixplot (ruwe of gestandaardiseerde data kan)
#ggplotMatrix(efa_data_filtered)


#Parallel analysis → aantal factoren bepalen
fa.parallel(efa_data_filtered, fa = "fa")

# EFA uitvoeren met gekozen aantal factoren (bv. 4)
fa_result <- fa(efa_data_filtered, nfactors = 25, rotate = "oblimin")
print(fa_result, cutoff = 0.3)
fa.diagram(fa_result, cut = 0.295, main = "Factor Analysis")

library(qgraph)
library(EGAnet)

# Gewone correlatiematrix
cor_matrix <- cor(efa_data_filtered, use = "pairwise.complete.obs")

# Plot netwerk gebaseerd op gewone correlaties
qgraph(cor_matrix, layout = "spring", labels = colnames(cor_matrix),
       title = "Pairwise correlation network", edge.color = "darkblue")

# Gebruik EGAnet om glasso netwerk op te bouwen
ega_network <- EGA(data = efa_data_filtered, model = "glasso", plot.EGA = FALSE)

# Extraheer netwerkmatrix
glasso_matrix <- ega_network$network

# Plot netwerk op basis van partiële correlaties
qgraph(glasso_matrix, layout = "spring", labels = colnames(glasso_matrix),
       title = "Regularized partial correlation network", edge.color = "darkred")

# 1. Count non-zero edges (off-diagonal)
nz_pairwise <- sum(cor_matrix[upper.tri(cor_matrix)] != 0)
nz_glasso   <- sum(glasso_matrix[upper.tri(glasso_matrix)] != 0)

message("Non-zero edges: pairwise = ", nz_pairwise,
        ", glasso = ", nz_glasso)

# 2. Edge-weight summaries
pw_edges  <- cor_matrix[upper.tri(cor_matrix)]
gl_edges  <- glasso_matrix[upper.tri(glasso_matrix)]

summary(pw_edges)
summary(gl_edges)

# 3. Centrality comparison
library(qgraph)
# compute strength centrality for each network
cent_pw    <- centrality_auto(qgraph(cor_matrix, DoNotPlot=TRUE))$node.centrality$Strength
cent_glass <- centrality_auto(qgraph(glasso_matrix, DoNotPlot=TRUE))$node.centrality$Strength

cor(cent_pw, cent_glass, method="spearman")

# 4. Adjacency similarity
# (Spearman between flattened upper triangles)
cor(pw_edges, gl_edges, method="spearman", use="pairwise.complete.obs")


#### PMFR (Pairwise Markov Random Field) #######
library(qgraph)
library(bootnet)
library(EGAnet)

# Estimate a Gaussian Graphical Model 
glasso_pmrf <- estimateNetwork(
  efa_data_filtered,
  default = "EBICglasso",
  corMethod = "cor_auto",      # Automatically chooses polychoric or Pearson
  tuning = 0.5                 # Regularization strength (default = 0.5)
)

# Plot the network
plot(glasso_pmrf,
     layout = "spring",
     title = "PMRF via EBICglasso",
     labels = colnames(efa_data_filtered))

# Centrality metrics (strength, betweenness, closeness)
centralityPlot(glasso_pmrf)

# Bootstrapping for accuracy/stability
pmrf_boot <- bootnet(
  glasso_pmrf,
  nBoots = 100,
  type = "case"
)

# Plot edge weight CIs
plot(pmrf_boot, labels = FALSE, order = "sample")

# Centrality stability (CS-coefficients)
cs_result <- corStability(pmrf_boot)
print(cs_result)

## Pmfr with Trustcountry measure ##
# Make a combined dataset
efa_pmrf_data <- efa_data_filtered %>%
  mutate(Trust_countrymeasure = data_clean$Trust_countrymeasure)

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



