"""This R code transforms the filtered data from RAP into long format (each row is a cognitive test). This prepares the data for GAMLSS"""


# Load required libraries
library(gamlss)
library(gamlss.dist)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

#import data
ukb_cog <- read_csv("ukb_cog_filtered.csv") #replace with data path

# Define cognitive tests
cognitive_tests <- c("PMat", "NM", "TMA", "TMB", "MP", "FI", "RT", "BL", "SD", "PAL", "TR", "PV")

# Simplified reshape function
reshape_for_gamlss <- function(ukb_cog) {
  
  test_data_list <- list()
  
  for(test in cognitive_tests) {
    cat("Processing test:", test, "\n")
    
    # Find all columns for this test
    test_cols <- grep(paste0("^", test, "-"), names(ukb_cog), value = TRUE)
    age_cols <- grep("^Age-", names(ukb_cog), value = TRUE)
    
    cat("Found columns for", test, ":", paste(test_cols, collapse = ", "), "\n")
    
    if(length(test_cols) > 0) {
      tryCatch({
        test_long <- ukb_cog %>%
          select(eid, Sex, Edu, all_of(test_cols), all_of(age_cols)) %>%
          pivot_longer(cols = all_of(test_cols), 
                       names_to = "test_instance", 
                       values_to = "test_score") %>%
          filter(!is.na(test_score)) %>%
          mutate(
            # Extract instance information
            instance = str_extract(test_instance, "(?<=-)[0-9]+$"),
            test_name = test,
            # Match age based on instance number
            age_at_test = case_when(
              instance == "0" ~ `Age-0`,
              instance == "1" ~ `Age-1`, 
              instance == "2" ~ `Age-2`,
              instance == "3" ~ `Age-3`,
              instance == "4" ~ `Age-4`,
              instance == "5" ~ `Age-5`,
              TRUE ~ NA_real_
            )
          ) %>%
          filter(!is.na(age_at_test)) %>%
          select(eid, Sex, Edu, test_name, instance, age_at_test, test_score)
        
        test_data_list[[test]] <- test_long
        cat("Created", nrow(test_long), "observations for", test, "\n")
        
        # Show instance breakdown
        if(nrow(test_long) > 0) {
          instance_summary <- test_long %>%
            count(instance) %>%
            arrange(instance)
          cat("Instances found:", paste(instance_summary$instance, collapse = ", "), "\n")
        }
        cat("\n")
        
      }, error = function(e) {
        cat("Error processing", test, ":", e$message, "\n\n")
      })
    } else {
      cat("No columns found for", test, "\n\n")
    }
  }
  
  # Combine all tests
  if(length(test_data_list) > 0) {
    combined_data <- bind_rows(test_data_list)
    cat("SUCCESS: Combined data created with", nrow(combined_data), "total observations\n")
    return(combined_data)
  } else {
    cat("No data was processed\n")
    return(NULL)
  }
}

# Run the reshape function
cat("Starting data reshape...\n")
long_data <- reshape_for_gamlss(ukb_cog)

# Remove NAs for Sex and Edu only (after reshape)
if(!is.null(long_data) && nrow(long_data) > 0) {
  
  cat("\n=== REMOVING NAs FOR SEX AND EDUCATION ===\n")
  cat("Before NA removal:", nrow(long_data), "observations\n")
  
  # Check how many NAs we have
  na_summary_before <- long_data %>%
    summarise(
      na_sex = sum(is.na(Sex)),
      na_edu = sum(is.na(Edu)),
      total_obs = n()
    )
  
  cat("NAs in Sex:", na_summary_before$na_sex, "\n")
  cat("NAs in Edu:", na_summary_before$na_edu, "\n")
  
  # Remove NAs for Sex and Edu only
  long_data_clean <- long_data %>%
    filter(
      !is.na(Sex),    # Remove participants missing sex
      !is.na(Edu)     # Remove participants missing education
    )
  
  cat("After NA removal:", nrow(long_data_clean), "observations\n")
  cat("Observations removed:", nrow(long_data) - nrow(long_data_clean), "\n")
  
  # Update long_data to the cleaned version
  long_data <- long_data_clean
  
  # Verify no NAs remain in Sex and Edu
  final_na_check <- long_data %>%
    summarise(
      na_sex = sum(is.na(Sex)),
      na_edu = sum(is.na(Edu))
    )
  
  if(final_na_check$na_sex == 0 && final_na_check$na_edu == 0) {
    cat("SUCCESS: No NAs remaining in Sex and Education\n")
  } else {
    cat("WARNING: Still have NAs after cleaning\n")
  }
}

# Verify the results
if(!is.null(long_data) && nrow(long_data) > 0) {
  cat("\n=== FINAL RESULTS ===\n")
  cat("Total observations:", nrow(long_data), "\n")
  cat("Unique participants:", length(unique(long_data$eid)), "\n")
  cat("Tests included:", paste(unique(long_data$test_name), collapse = ", "), "\n")
  
  # Show instance distribution
  cat("\n=== INSTANCE DISTRIBUTION ===\n")
  instance_dist <- long_data %>%
    count(test_name, instance) %>%
    spread(instance, n, fill = 0)
  print(instance_dist)
  write.csv(instance_dist, "instance_distribution.csv", row.names = FALSE)
  
  # Show age ranges
  cat("\n=== AGE SUMMARY ===\n")
  age_summary <- long_data %>%
    group_by(instance) %>%
    summarise(
      min_age = min(age_at_test, na.rm = TRUE),
      max_age = max(age_at_test, na.rm = TRUE),
      mean_age = round(mean(age_at_test, na.rm = TRUE), 1),
      .groups = 'drop'
    )
  print(age_summary)
  write.csv(age_summary, "age_summary.csv", row.names = FALSE)
  
  # Sample of the data
  cat("\n=== SAMPLE DATA ===\n")
  print(head(long_data, 10))
  
  # Save for safety
  save(long_data, file = "long_data.RData")
  cat("\nData saved to long_data.RData\n")
  
} else {
  cat("FAILED: long_data was not created or is empty\n")
}



