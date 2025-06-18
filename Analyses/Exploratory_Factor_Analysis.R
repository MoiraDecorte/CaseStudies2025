# Loading packages
library(readr)
library(dplyr)
library(psych)
library(lme4)
library(lmerTest)

# reading and loading in the data
raw_data <- read_csv("C:/Users/moira/Downloads/COVIDiSTRESS global survey May 30 2020 (___final cleaned file___)/COVIDiSTRESS global survey May 30 2020 (___final cleaned file___).csv",
                     locale = locale(encoding = "ISO-8859-1"), na = c("NA", ""))

names(raw_data) <- gsub("[ -]", "_", names(raw_data))

# Selecting all items using a Likert-scale 
likert_items <- grep("^(Scale_PSS10_UCLA_|Corona_concerns_|Compliance_|BFF_15_|Expl_Distress_|SPS_|Expl_Coping_|Expl_media_)", names(raw_data), value = TRUE)

# Deleting "Expl_Distress_txt" as it is a textfield, and cannot be part of an EFA
likert_items <- likert_items[likert_items != "Expl_Distress_txt"]

# Separating the SPS group due to low correlations with other items, which could distort the EFA results
sps_items <- grep("^SPS_", likert_items, value = TRUE)

# Delete SPS-items from EFA data set 
likert_items_no_sps <- likert_items[!(likert_items %in% sps_items)]

# Grouping items 
oecd_items <- grep("^(OECD_people_|OECD_insititutions_)", names(raw_data), value = TRUE)
extra_items <- c("Trust_countrymeasure")
all_items <- c(likert_items_no_sps, oecd_items, extra_items)

# Cleaning the data before EFA 
efa_data <- raw_data[, all_items]
efa_data[efa_data == 99] <- NA # Replace '99' values (used as missing data placeholder) with NA
efa_data[] <- lapply(efa_data, function(x) as.numeric(as.character(x)))

# 7 means "Does not apply" in the Expl_Distress-questions so we put it at NA 
expl_distress_items <- grep("^Expl_Distress_", names(efa_data), value = TRUE)
efa_data[, expl_distress_items][efa_data[, expl_distress_items] == 7] <- NA

# Add country and subjectID 
efa_data$country <- raw_data$Country
efa_data$SubjectID <- as.factor(1:nrow(efa_data))

# Standardizing all scales besides Trust_countrymeasure and removes rows with missing country variable 
data_scaled <- efa_data %>%
  mutate(across(where(is.numeric) & !any_of(c("Trust_countrymeasure")), ~ scale(.))) %>%
  filter(!is.na(country))

# Keeping only rows with at least 5 non-missing items
data_efa <- data_scaled %>%
  filter(rowSums(!is.na(select(., -country, -Trust_countrymeasure, -SubjectID))) >= 5)

# preparing data for EFA by only selecting the relevant variables
fa_input <- select(data_efa, -country, -Trust_countrymeasure, -SubjectID)

# Correlationmatrix
cor_matrix <- cor(fa_input, use = "pairwise.complete.obs")

# Items with max correlatie < .30
low_correlation_items <- names(which(apply(abs(cor_matrix), 2, function(x) max(x[-which.max(x)], na.rm = TRUE)) < 0.3))

# Print these items
print("Items met te lage correlatie:")
print(low_correlation_items)

# Remove them from dataset
fa_input_final <- fa_input[, !(names(fa_input) %in% low_correlation_items)]

# Calculate the KMO
final_kmo <- KMO(fa_input_final)
print(final_kmo)

# Bartlett on complete cases
fa_input_complete_final <- fa_input_final[complete.cases(fa_input_final), ]
bartlett_final <- cortest.bartlett(fa_input_complete_final)
print(bartlett_final)

# executing EFA — Scree plot + Parallel Analysis + EFA
fa.parallel(fa_input_final, fa = "fa", n.iter = 20, show.legend = TRUE, main = "Parallel Analysis Scree Plot")
n_factors <- 9 # based on the elbow/knee in screeplot
fa_result <- fa(fa_input_final, nfactors = n_factors, rotate = "oblimin", scores = "regression", fm = "minres")
print(fa_result)

# Visualising factor structure
fa.diagram(fa_result)

# save as pdf to be able to see all the variables:
pdf("efa_factor_diagram.pdf", width = 35, height = 20)
fa.diagram(fa_result, cut = 0.3)  # cut minimum loading being shown
dev.off()

# Load pattern matrix in a data.frame
loadings_matrix <- as.data.frame(fa_result$loadings[1:nrow(fa_result$loadings), ])
loadings_matrix$Communality <- fa_result$communality

# Vector with factors in ascending order  
desired_order <- paste0("MR", 1:n_factors) # 'MR' = Manually labeled factor numbers (e.g., MR1 = first factor based on EFA loadings)
existing_factors <- intersect(desired_order, colnames(loadings_matrix))
loadings_ordered <- loadings_matrix[, c(existing_factors, "Communality")]

# Round with 3 decimals
library(dplyr)

loadings_ordered <- loadings_ordered %>%
  mutate(across(all_of(existing_factors), ~ round(., 3)),
         Communality = round(Communality, 3))
head(loadings_ordered)

# saving as CSV to have the table for our paper 
write.csv(loadings_ordered, "efa_factor_loadings_ordered_2.csv")

# Summary of the explained variance per factor 
explained_variance <- data.frame(
  Factor = paste0("MR", 1:n_factors),
  SS_Loadings = fa_result$Vaccounted["SS loadings", ],
  Proportion_Variance = fa_result$Vaccounted["Proportion Var", ],
  Cumulative_Variance = fa_result$Vaccounted["Cumulative Var", ]
)
print(explained_variance)

# Minimum and maximum communality
min_comm <- min(fa_result$communality, na.rm = TRUE)
max_comm <- max(fa_result$communality, na.rm = TRUE)
cat("Minimum communality:", round(min_comm, 3), "\n")
cat("Maximum communality:", round(max_comm, 3), "\n")

MR1_items <- c("Compliance_4", 
               paste0("Expl_Distress_", c(1:9, 11, 15:20, 22:24)))

MR2_items <- c(paste0("BFF_15_", 4:6),
               paste0("Expl_Coping_", c(2:5, 10)))

MR3_items <- c("Expl_Coping_1", "Expl_media_1", 
               "OECD_people_1", "OECD_people_2", 
               paste0("OECD_insititutions_", 1:6))

MR4_items <- c(paste0("Corona_concerns_", 1:5), 
               paste0("Expl_Distress_", 12:14))

MR5_items <- c(paste0("Expl_Coping_", c(10:12, 15, 16)))

MR6_items <- c(paste0("Scale_PSS10_UCLA_", 1:10), 
               paste0("BFF_15_", c(1:3, 7:9, 12:15)))

MR7_items <- c("Expl_Coping_7", "Expl_Coping_8")

MR8_items <- c(paste0("Expl_media_", 1:5))

MR9_items <- c(paste0("Compliance_", 1:3), "BFF_15_3")

alpha_results <- list()

for (i in 1:9) {
  item_vector_name <- paste0("MR", i, "_items")
  
  if (exists(item_vector_name)) {
    items <- get(item_vector_name)
    data_subset <- fa_input_final[, items]
    
    # Calculate alpha with auto-reverse check
    alpha_result <- psych::alpha(data_subset, check.keys = TRUE)
    
    # Store full result
    alpha_results[[paste0("MR", i)]] <- alpha_result
    
    # Print concise summary
    cat("\n====================\n")
    cat("Factor:", paste0("MR", i), "\n")
    print(alpha_result$total[["raw_alpha"]])
    cat("====================\n")
  } else {
    warning(paste("Item vector", item_vector_name, "not found."))
  }
}

factor_scores <- as.data.frame(fa_result$scores)

# Create table with eigenvalues, % variance, and cumulative %
explained_variance_table <- data.frame(
  Factor = paste0("MR", 1:n_factors),
  Eigenvalue = round(fa_result$Vaccounted["SS loadings", ], 2),
  Percent_of_Variance = round(fa_result$Vaccounted["Proportion Var", ] * 100, 2),
  Cumulative_Percent = round(fa_result$Vaccounted["Cumulative Var", ] * 100, 2)
)

# Print table to console
print(explained_variance_table)                                           

# Add Trust, country and SubjectID for later analysis
factor_scores$Trust_countrymeasure <- data_efa$Trust_countrymeasure
factor_scores$country <- data_efa$country
factor_scores$SubjectID <- data_efa$SubjectID

# Add continent via mapping
factor_scores$country <- tolower(trimws(factor_scores$country))

continent_lookup <- bind_rows(
  data.frame(country = tolower(c("The Gambia", "Tanzania", "South Sudan", "Sudan, South", "Guinea", "Guinea-Bissau", "Ghana", "Gabon", "Ethiopia", "Eritrea", "Egypt", "Djibouti", "Côte d’Ivoire", "Congo, Republic of the", "Congo, Democratic Republic of the", "Comoros", "Cameroon", "Cabo Verde", "Burkina Faso", "Burundi", "Benin", "Angola", "Algeria", "Kenya", "Lesotho", "Libya", "Madagascar", "Malawi", "Mali", "Mauritania", "Mauritius", "Morocco", "Mozambique", "Namibia", "Niger", "Nigeria", "Rwanda", "Sao Tome and Principe", "Senegal", "Seychelles", "Somalia", "South Africa", "Sudan", "Tunisia", "Uganda", "Zambia", "Zimbabwe")), Continent = "Africa"),
  data.frame(country = tolower(c("East Timor (Timor-Leste)", "China", "India", "Japan", "South Korea", "North Korea", "Nepal", "Bangladesh", "Pakistan", "Thailand", "Vietnam", "Indonesia", "Philippines", "Malaysia", "Singapore", "Sri Lanka", "Afghanistan", "Iran", "Iraq", "Jordan", "Kazakhstan", "Kuwait", "Lebanon", "Myanmar (Burma)", "Oman", "Qatar", "Saudi Arabia", "Syria", "Taiwan", "United Arab Emirates", "Uzbekistan")), Continent = "Asia"),
  data.frame(country = tolower(c("Austria", "Belgium", "Croatia", "Denmark", "Finland", "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Italy", "Netherlands", "Norway", "Poland", "Portugal", "Spain", "Sweden", "Switzerland", "United Kingdom")), Continent = "Europe"),
  data.frame(country = tolower(c("Canada", "United States", "Mexico", "Cuba", "Dominican Republic", "Guatemala", "Honduras", "Panama")), Continent = "North America"),
  data.frame(country = tolower(c("Argentina", "Brazil", "Chile", "Colombia", "Peru", "Uruguay", "Venezuela")), Continent = "South America"),
  data.frame(country = tolower(c("Australia", "New Zealand", "Fiji", "Papua New Guinea")), Continent = "Oceania")
)

factor_scores <- left_join(factor_scores, continent_lookup, by = "country")

# Save
saveRDS(factor_scores, "factor_scores.rds")
saveRDS(fa_result, "fa_result.rds")
saveRDS(data_efa, "data_efa.rds")
