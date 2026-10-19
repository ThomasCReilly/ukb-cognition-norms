"""In preparation for normative modelling, this code was used to calculate the data shape and structure and present it as a clean table for viewing.
It will also create a series of graphs to visualise the raw data. """

# Load required libraries
library(ggplot2)
library(gridExtra)
library(dplyr)

# Complete comprehensive data summary for distribution selection
# May need to import and define long_data if not already defined from previous script

create_distribution_summary <- function(long_data) {
  
  # Calculate skewness manually
  calc_skewness <- function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    if(n < 3) return(NA)
    m3 <- mean((x - mean(x))^3)
    m2 <- mean((x - mean(x))^2)
    skewness <- m3 / (m2^(3/2))
    return(skewness)
  }
  
  # Calculate kurtosis manually
  calc_kurtosis <- function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    if(n < 4) return(NA)
    m4 <- mean((x - mean(x))^4)
    m2 <- mean((x - mean(x))^2)
    kurtosis <- m4 / (m2^2)
    return(kurtosis)
  }
  
  # Calculate outliers using IQR method - separate above and below
  detect_outliers_separate <- function(x) {
    Q1 <- quantile(x, 0.25, na.rm = TRUE)
    Q3 <- quantile(x, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    
    below_outliers <- sum(x < lower_bound, na.rm = TRUE)
    above_outliers <- sum(x > upper_bound, na.rm = TRUE)
    total_outliers <- below_outliers + above_outliers
    
    return(list(
      below = below_outliers,
      above = above_outliers,
      total = total_outliers
    ))
  }
  
  # Normality test (sample for large datasets)
  safe_shapiro <- function(x) {
    x_clean <- x[!is.na(x)]
    if(length(x_clean) < 3) return(NA)
    if(length(x_clean) > 5000) x_clean <- sample(x_clean, 5000)
    tryCatch(shapiro.test(x_clean)$p.value, error = function(e) NA)
  }
  
  summary_table <- long_data %>%
    group_by(test_name) %>%
    summarise(
      n_obs = n(),
      mean = round(mean(test_score, na.rm = TRUE), 3),
      median = round(median(test_score, na.rm = TRUE), 3),
      Q1 = round(quantile(test_score, 0.25, na.rm = TRUE), 3),
      Q3 = round(quantile(test_score, 0.75, na.rm = TRUE), 3),
      sd = round(sd(test_score, na.rm = TRUE), 3),
      min_val = round(min(test_score, na.rm = TRUE), 3),
      max_val = round(max(test_score, na.rm = TRUE), 3),
      range_span = round(max_val - min_val, 3),
      skewness = round(calc_skewness(test_score), 3),
      kurtosis = round(calc_kurtosis(test_score), 3),
      excess_kurtosis = round(calc_kurtosis(test_score) - 3, 3),
      prop_zeros = round(mean(test_score == 0, na.rm = TRUE), 4),
      n_zeros = sum(test_score == 0, na.rm = TRUE),
      prop_negative = round(mean(test_score < 0, na.rm = TRUE), 4),
      n_negative = sum(test_score < 0, na.rm = TRUE),
      outliers = list(detect_outliers_separate(test_score)),
      shapiro_p = safe_shapiro(test_score),
      is_normal_05 = ifelse(is.na(shapiro_p), "Unknown", 
                            ifelse(shapiro_p > 0.05, "Yes", "No")),
      is_integer = all(test_score == round(test_score, 0), na.rm = TRUE),
      is_bounded_01 = all(test_score >= 0 & test_score <= 1, na.rm = TRUE),
      is_positive = all(test_score >= 0, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      # Extract outlier counts from list column
      n_outliers_below = sapply(outliers, function(x) x$below),
      n_outliers_above = sapply(outliers, function(x) x$above),
      n_outliers_total = sapply(outliers, function(x) x$total),
      prop_outliers = round(n_outliers_total / n_obs, 4),
      skew_interp = case_when(
        is.na(skewness) ~ "Unknown",
        abs(skewness) < 0.5 ~ "Symmetric",
        skewness >= 0.5 & skewness < 1 ~ "Moderate Right",
        skewness >= 1 ~ "High Right", 
        skewness <= -0.5 & skewness > -1 ~ "Moderate Left",
        skewness <= -1 ~ "High Left"
      ),
      kurt_interp = case_when(
        is.na(excess_kurtosis) ~ "Unknown",
        abs(excess_kurtosis) < 0.5 ~ "Normal tails",
        excess_kurtosis >= 0.5 ~ "Heavy tails",
        excess_kurtosis <= -0.5 ~ "Light tails"
      ),
      recommended_families = case_when(
        prop_zeros > 0.1 ~ "Zero-inflated (ZAGA, ZAIG, ZANBI)",
        is_bounded_01 ~ "Beta (BE, BEINF)",
        !is_positive ~ "Normal-based (NO, TF, ST1, SN1)",
        is.na(skewness) | is.na(excess_kurtosis) ~ "Normal-based (NO, TF)",
        abs(skewness) < 0.5 & abs(excess_kurtosis) < 0.5 ~ "Normal-based (NO, TF)",
        abs(skewness) >= 0.5 & abs(excess_kurtosis) < 0.5 ~ "Skewed (SN1, ST1)",
        abs(skewness) < 0.5 & abs(excess_kurtosis) >= 0.5 ~ "Heavy-tailed (TF, ST1)",
        TRUE ~ "Flexible (ST1, ST2, BCCG)"
      )
    ) %>%
    # Remove the list column
    select(-outliers)
  
  return(summary_table)
}

# Create the clean table
create_clean_distribution_table <- function(distribution_summary) {
  clean_table <- distribution_summary %>%
    mutate(
      N = formatC(n_obs, format="d", big.mark=","),
      `Mean±SD` = paste0(sprintf("%.2f", mean), "±", sprintf("%.2f", sd)),
      Median = sprintf("%.2f", median),
      `Q1-Q3` = paste0(sprintf("%.2f", Q1), " - ", sprintf("%.2f", Q3)),
      Range = paste0(min_val, " to ", max_val),
      `Skew (Type)` = paste0(sprintf("%.2f", skewness), " (", skew_interp, ")"),
      `Kurt (Type)` = paste0(sprintf("%.2f", excess_kurtosis), " (", kurt_interp, ")"),
      Zeros = paste0(n_zeros, " (", sprintf("%.1f", prop_zeros*100), "%)"),
      Negs = ifelse(n_negative > 0, paste0(n_negative, " (", sprintf("%.1f", prop_negative*100), "%)"), "None"),
      `Out Below` = as.character(n_outliers_below),
      `Out Above` = as.character(n_outliers_above),
      `Out Total` = paste0(n_outliers_total, " (", sprintf("%.1f", prop_outliers*100), "%)")
    ) %>%
    select(
      Test = test_name,
      N,
      `Mean±SD`,
      Median,
      `Q1-Q3`,
      Range,
      `Skew (Type)`,
      `Kurt (Type)`, 
      Zeros,
      Negs,
      `Out Below`,
      `Out Above`,
      `Out Total`,
      Normal = is_normal_05,
      `Recommended Families` = recommended_families
    )
  
  return(clean_table)
}

# Function to detect outliers using IQR method for plotting
detect_outliers_for_plots <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  list(
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    outlier_indices = which(x < lower_bound | x > upper_bound)
  )
}

# Function to create comprehensive plots for each test
plot_cognitive_tests <- function(long_data) {
  
  # Get unique tests
  tests <- unique(long_data$test_name)
  
  # Create main directory
  dir.create("cogtest_descriptive_analysis", showWarnings = FALSE)
  
  all_outliers <- list()
  
  for(test in tests) {
    
    test_data <- long_data %>% filter(test_name == test)
    
    # Create test directory
    test_dir <- file.path("cogtest_descriptive_analysis", test)
    dir.create(test_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Detect outliers
    outlier_info <- detect_outliers_for_plots(test_data$test_score)
    
    # Create outlier datasets with all columns
    below_outliers <- test_data[test_data$test_score < outlier_info$lower_bound, ]
    above_outliers <- test_data[test_data$test_score > outlier_info$upper_bound, ]
    
    # Store outliers for this test
    all_outliers[[test]] <- list(
      below_range = below_outliers,
      above_range = above_outliers,
      lower_bound = outlier_info$lower_bound,
      upper_bound = outlier_info$upper_bound
    )
    
    # Add outlier flag to test data
    test_data$is_outlier <- test_data$test_score < outlier_info$lower_bound | 
      test_data$test_score > outlier_info$upper_bound
    
    # Plot 1: Histogram with density curve
    cat("  Creating histogram...")
    p1 <- ggplot(test_data, aes(x = test_score)) +
      geom_histogram(aes(y = after_stat(density)), bins = 50, 
                     fill = "grey90", color = "grey70", alpha = 0.8, linewidth = 0.05) +
      geom_density(color = "grey40", linewidth = 0.15) +
      labs(title = paste(test, "- Distribution"), 
           x = "Test Score", y = "Density") +
      theme_minimal()
    
    ggsave(file.path(test_dir, paste0(test, "_histogram.png")), p1, 
           width = 6, height = 4, dpi = 600, bg = "white")
    cat(" saved\n")
    
    # Plot 2: Q-Q plot for normality assessment
    cat("  Creating Q-Q plot...")
    p2 <- ggplot(test_data, aes(sample = test_score)) +
      stat_qq(color = "grey70", size = 0.15, alpha = 0.5) +
      stat_qq_line(color = "grey30", linetype = "dashed", linewidth = 0.2) +
      labs(title = paste(test, "- Q-Q Plot"), 
           x = "Theoretical Quantiles", y = "Sample Quantiles") +
      theme_minimal()
    
    ggsave(file.path(test_dir, paste0(test, "_qqplot.png")), p2, 
           width = 6, height = 4, dpi = 600, bg = "white")
    cat(" saved\n")
    
    # Plot 3: Horizontal box plot with outliers highlighted
    cat("  Creating box plot...")
    p3 <- ggplot(test_data, aes(x = test_score, y = "")) +
      geom_boxplot(fill = "grey90", color = "grey10", alpha = 0.8, 
                   width = 0.4, linewidth = 0.25, outlier.color = "grey20", 
                   outlier.size = 0.4, outlier.alpha = 0.5) +
      geom_point(data = test_data[test_data$is_outlier, ], 
                 aes(x = test_score, y = ""), 
                 color = "grey20", size = 0.4, alpha = 0.6, shape = 1) +
      labs(title = paste(test, "- Box Plot"), 
           x = "Test Score", y = "") +
      theme_minimal() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    
    ggsave(file.path(test_dir, paste0(test, "_boxplot.png")), p3, 
           width = 8, height = 3, dpi = 600, bg = "white")
    cat(" saved\n")
    
    # Plot 4: Violin plot with data points and outliers
    p4 <- ggplot(test_data, aes(x = test_score, y = "")) +
      geom_jitter(aes(shape = is_outlier), width = 0.15, size = 0.15, 
                  alpha = 0.3, color = "grey80") +
      geom_violin(fill = "grey90", color = "grey10", alpha = 0.8, linewidth = 0.15) +
      scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 1), 
                         labels = c("Normal", "Outlier"),
                         name = "Points") +
      labs(title = paste(test, "- Violin Plot"), 
           x = "Test Score", y = "") +
      theme_minimal()
    
    ggsave(file.path(test_dir, paste0(test, "_violin.png")), p4, 
           width = 8, height = 5, dpi = 600, bg = "white")
    
    # Plot 5: Score by age (scatter with smooth)
    p5 <- ggplot(test_data, aes(x = test_score, y = age_at_test)) +
      geom_point(aes(shape = is_outlier), alpha = 0.25, size = 0.35, color = "grey40") +
      scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16), 
                         labels = c("Normal", "Outlier"),
                         name = "Points") +
      labs(title = paste(test, "- Score vs Age"), 
           x = "Test Score", y = "Age") +
      theme_minimal()
    
    ggsave(file.path(test_dir, paste0(test, "_age_scatter.png")), p5, 
           width = 10, height = 5, dpi = 600, bg = "white")
    
    # Print outlier summary for this test
    cat("\n", rep("=", 60), "\n")
    cat("OUTLIER SUMMARY:", test, "\n")
    cat(rep("=", 60), "\n")
    cat("Bounds: [", round(outlier_info$lower_bound, 3), ", ", 
        round(outlier_info$upper_bound, 3), "]\n")
    cat("Below range:", nrow(below_outliers), "| Above range:", nrow(above_outliers), "\n")
    cat("Plots saved to:", test_dir, "\n")
    
    if(nrow(below_outliers) > 0) {
      cat("\nBELOW RANGE OUTLIERS:\n")
      print(below_outliers[, c("eid", "test_score", "age_at_test", "Sex", "Edu")])
    }
    
    if(nrow(above_outliers) > 0) {
      cat("\nABOVE RANGE OUTLIERS:\n")
      print(above_outliers[, c("eid", "test_score", "age_at_test", "Sex", "Edu")])
    }
    
    # Save outlier data to CSV files in test directory
    if(nrow(below_outliers) > 0) {
      write.csv(below_outliers, 
                file.path(test_dir, paste0(test, "_below_outliers.csv")), 
                row.names = FALSE)
    }
    
    if(nrow(above_outliers) > 0) {
      write.csv(above_outliers, 
                file.path(test_dir, paste0(test, "_above_outliers.csv")), 
                row.names = FALSE)
    }
  }
  
  cat("\n", rep("=", 60), "\n")
  cat("ALL PLOTS COMPLETED\n")
  cat("Directory structure: cogtest_descriptive_analysis/[test_name]/\n")
  cat("Plot files: histogram, qqplot, boxplot, violin, age_scatter\n")
  cat("Outlier files: [test]_below_outliers.csv, [test]_above_outliers.csv\n")
  cat(rep("=", 60), "\n")
  
  return(all_outliers)
}

# MAIN EXECUTION
cat("=== COGNITIVE TEST DISTRIBUTION ANALYSIS ===\n\n")

# Generate the distribution summary first
cat("1. Generating distribution summary...\n")

# Create main directory first
dir.create("cogtest_descriptive_analysis", showWarnings = FALSE)

distribution_summary <- create_distribution_summary(long_data)
clean_dist_table <- create_clean_distribution_table(distribution_summary)

# Print with nice formatting
cat("\n")
cat(rep("=", 120), "\n")
cat("COGNITIVE TEST DISTRIBUTION SUMMARY FOR GAMLSS FAMILY SELECTION\n")
cat(rep("=", 120), "\n\n")

# Print each row with formatting
for(i in 1:nrow(clean_dist_table)) {
  row <- clean_dist_table[i,]
  
  cat("TEST:", sprintf("%-6s", row$Test), "| N:", sprintf("%8s", row$N), 
      "| Mean±SD:", sprintf("%12s", row$`Mean±SD`), 
      "| Median:", sprintf("%8s", row$Median), 
      "| Q1-Q3:", sprintf("%15s", row$`Q1-Q3`), "\n")
  
  cat("     Shape: Skew:", sprintf("%-20s", row$`Skew (Type)`), 
      "| Kurt:", sprintf("%-18s", row$`Kurt (Type)`), "\n")
  
  cat("     Issues: Zeros:", sprintf("%-15s", row$Zeros), 
      "| Negs:", sprintf("%-15s", row$Negs), 
      "| Out Below:", sprintf("%-8s", row$`Out Below`), 
      "| Out Above:", sprintf("%-8s", row$`Out Above`), "\n")
  
  cat("     Tests: Normal:", sprintf("%-8s", row$Normal), "\n")
  
  cat("     RECOMMENDED:", row$`Recommended Families`, "\n")
  
  cat(rep("-", 120), "\n")
}

# Also print as a regular table
cat("\n\nTABULAR FORMAT:\n")
cat(rep("=", 120), "\n")
print(clean_dist_table)

# Save distribution summary to CSV
write.csv(clean_dist_table, "cogtest_descriptive_analysis/cognitive_test_distribution_summary.csv", row.names = FALSE)
cat("\nDistribution summary saved to: cogtest_descriptive_analysis/cognitive_test_distribution_summary.csv\n")

# Generate plots
cat("\n2. Generating diagnostic plots...\n")
all_outliers <- plot_cognitive_tests(long_data)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files created:\n")
cat("- cogtest_descriptive_analysis/cognitive_test_distribution_summary.csv (summary table)\n")
cat("- cogtest_descriptive_analysis/[test_name]/ (individual test directories with plots and outlier data)\n")