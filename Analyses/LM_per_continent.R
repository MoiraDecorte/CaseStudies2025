# Loading packages
library(readr)
library(dplyr)
library(psych)
library(lme4)
library(lmerTest)

# Loading datasets
factor_scores <- readRDS("factor_scores.rds")
data_efa <- readRDS("data_efa.rds")
fa_result <- readRDS("fa_result.rds")

# LM per continent
model_results <- list()
continents <- sort(unique(factor_scores$Continent))

for (cont in continents) {
  cat("\n===============================\n")
  cat("\U0001F4CD Resultaten voor continent:", cont, "\n")
  cat("===============================\n")
  flush.console()
  
  data_sub <- factor_scores %>% filter(Continent == cont)
  
  complete_rows <- complete.cases(select(data_sub, starts_with("MR")), data_sub$Trust_countrymeasure)
  data_sub <- data_sub[complete_rows, ]
  
  if (nrow(data_sub) < 10) {
    cat("\u26A0\uFE0F Te weinig bruikbare data voor ", cont, "(", nrow(data_sub), "rijen)\n")
    next
  }
  
  predictors <- paste0(colnames(data_sub)[grepl("^MR", colnames(data_sub))], collapse = " + ")
  formula <- as.formula(paste("Trust_countrymeasure ~", predictors, "+ (1 | SubjectID)"))
  
  tryCatch({
    model <- lmer(formula, data = data_sub)
    model_results[[cont]] <- model
    print(summary(model))
  }, error = function(e) {
    cat("\u2139\uFE0F Enkelvoudige observaties per subject \u2192 gebruik lm()\n")
    formula_simple <- as.formula(paste("Trust_countrymeasure ~", predictors))
    model <- lm(formula_simple, data = data_sub)
    model_results[[cont]] <- model
    print(summary(model))
  })
}