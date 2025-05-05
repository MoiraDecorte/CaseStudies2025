library(readr)

# Read the csv file
data <- read_csv(
  "C:\\Users\\piepe\\OneDrive - UGent\\UGENT\\2024 - 2025\\CASE STUDIES\\DATA\\COVIDiSTRESS global survey May 30 2020 (___final cleaned file___).csv",
  locale = locale(encoding = "UTF-8"),
  na = c("NA", "")
)

# Eerste paar rijen bekijken
head(data)
colnames(data)