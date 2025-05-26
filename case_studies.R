#Packages laden
library(readr)
library(dplyr)
library(psych)
library(GGally)
library(ggplot2)
library(gridExtra)

#Data inlezen
data_raw <- read_csv(
  "C:/Users/moira/Downloads/COVIDiSTRESS global survey May 30 2020 (___final cleaned file___)/COVIDiSTRESS global survey May 30 2020 (___final cleaned file___).csv",
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

efa_data <- data_clean[, likert_items]

#'7 = niet van toepassing' omzetten naar NA
efa_data[efa_data == 7] <- NA

#Verwijder variabelen zonder voldoende correlatie
remove_uncorrelated_vars <- function(df, threshold = 0.3) {
  numeric_cols <- sapply(df, is.numeric)
  df <- df[, numeric_cols]
  cor_matrix <- cor(df, use = "pairwise.complete.obs")
  has_correlations <- sapply(1:ncol(cor_matrix), function(i) {
    any(abs(cor_matrix[i, -i]) >= threshold)
  })
  df_filtered <- df[, has_correlations]
  return(df_filtered)
}

efa_data_filtered <- remove_uncorrelated_vars(efa_data, threshold = 0.3)

#Visualiseer correlaties
ggplotMatrix <- function(efa_data, bin_width = 0.5, axis_range = c(1, 6)) {
  my_bins <- seq(axis_range[1], axis_range[2], bin_width)
  axis_breaks <- axis_range[1]:axis_range[2]
  
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


cor_matrix <- cor(efa_data_filtered, use = "pairwise.complete.obs")
heatmap(cor_matrix)

library(corrplot)
corrplot(cor_matrix, method = "color", tl.cex = 0.5)
removed_vars <- setdiff(names(efa_data), names(efa_data_filtered))
removed_vars
efa_numeric <- efa_data[, sapply(efa_data, is.numeric)]
cor_matrix <- cor(efa_numeric, use = "pairwise.complete.obs")

removed_vars_numeric <- intersect(removed_vars, colnames(efa_numeric))

cor_details <- data.frame(
  variable = removed_vars_numeric,
  max_correlation = sapply(removed_vars_numeric, function(var) {
    max(abs(cor_matrix[var, setdiff(colnames(cor_matrix), var)]), na.rm = TRUE)
  })
)

# Sorteer oplopend om te zien welke écht zwak waren
cor_details <- cor_details[order(cor_details$max_correlation), ]
print(cor_details)



#Parallel analysis → aantal factoren bepalen
fa.parallel(efa_data_filtered, fa = "fa")

# EFA uitvoeren met gekozen aantal factoren (bv. 4)
fa_result <- fa(efa_data_filtered, nfactors = 20, rotate = "oblimin")
print(fa_result, cutoff = 0.3)
fa.diagram(fa_result, cut = 0.295, main = "Factor Analysis")


library(EGAnet)
ega_result <- EGA(data = efa_data_filtered, model = "glasso", plot.EGA = TRUE)
efa_data_filtered[efa_data_filtered == 99] <- NA

