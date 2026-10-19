#This script transforms the Prospective Memory data into long format to be re-integrated into the dataset"

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# Load the PMem data
PMem_data <- read_csv("PMem_data.csv")

# Transform PMem to long format
PMem_long <- PMem_data %>%
  select(eid, Sex, Edu, starts_with("Age-"), starts_with("PMem-")) %>%
  pivot_longer(
    cols = starts_with("PMem-"),
    names_to = "test_instance",
    values_to = "test_score"
  ) %>%
  filter(!is.na(test_score)) %>%
  mutate(
    # Extract instance number
    instance = str_extract(test_instance, "(?<=-)[0-9]+$"),
    test_name = "PMem",
    # Match age based on instance number
    age_at_test = case_when(
      instance == "0" ~ `Age-0`,
      instance == "1" ~ `Age-1`,
      instance == "2" ~ `Age-2`,
      instance == "3" ~ `Age-3`,
      instance == "4" ~ `Age-4`,
      instance == "5" ~ `Age-5`,
      TRUE ~ NA_real_
    ),
    # PMem is already 0/1, so use as impaired directly
    impaired = test_score,
    # Set z_score and percentile to NA (since PMem is binary, not continuous)
    z_score = NA_real_,
    percentile = NA_real_
  ) %>%
  filter(!is.na(age_at_test)) %>%
  # Remove NAs for Sex and Edu to match existing data
  filter(!is.na(Sex), !is.na(Edu)) %>%
  # Select columns in the same order as cv_results_all_test_cleaned.csv
  select(eid, Sex, Edu, test_name, instance, age_at_test, test_score, z_score, percentile, impaired)

# Check the PMem data
cat("\n=== PMem DATA SUMMARY ===\n")
cat("PMem observations:", nrow(PMem_long), "\n")
cat("Unique participants:", length(unique(PMem_long$eid)), "\n")

# Show instance breakdown
cat("\n=== INSTANCE DISTRIBUTION ===\n")
instance_summary <- PMem_long %>%
  count(instance) %>%
  arrange(instance)
print(instance_summary)

# Show impairment distribution (0s and 1s)
cat("\n=== IMPAIRMENT DISTRIBUTION ===\n")
impaired_dist <- PMem_long %>%
  count(impaired) %>%
  mutate(percentage = round(100 * n / sum(n), 2))
print(impaired_dist)

# Show age ranges by instance
cat("\n=== AGE SUMMARY BY INSTANCE ===\n")
age_summary <- PMem_long %>%
  group_by(instance) %>%
  summarise(
    min_age = min(age_at_test, na.rm = TRUE),
    max_age = max(age_at_test, na.rm = TRUE),
    mean_age = round(mean(age_at_test, na.rm = TRUE), 1),
    n = n(),
    .groups = 'drop'
  )
print(age_summary)

# Sample of the data
cat("\n=== SAMPLE DATA ===\n")
print(head(PMem_long, 10))

# Save the PMem_long data
save(PMem_long, file = "PMem_long.RData")
write_csv(PMem_long, "PMem_long.csv")
cat("\nPMem long data saved to PMem_long.RData and PMem_long.csv\n")