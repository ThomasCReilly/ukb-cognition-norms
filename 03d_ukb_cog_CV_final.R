
# GAMLSS Cross-Validation for UK Biobank Cognitive Tests

# THis is the script that fits the optimal distirbution and uses a cross validation design
# to generate z-scores. For each cognitive test, it fits the optimal GAMLSS distribution (selected by
# prior dist_test) using 10-fold cross-validation by participant. It uses only baseline data for training
# but then uses repeated measures for prediction. Generates
# out-of-sample z-scores and percentiles for downstream analyses. 
#
# MODEL STRUCTURE
#   mu     ~ pb(age) + Sex + pb(Edu)
#   sigma  ~ pb(age) + Sex + pb(Edu)
#   nu     ~ 1
#   tau    ~ 1
#
#   The random() term is a participant-level random intercept that accounts for
#   within-participant correlation across repeated cognitive assessments.
#
# CV STRATEGY
#   10-fold cross-validation, splitting on participant ID (eid) so all visits
#   for a given participant fall in the same fold. Predicted distribution
#   parameters are used with the CDF (p*) to compute test-set percentiles,
#   then converted to z-scores via qnorm().
#
# DATA REQUIREMENTS
#   long_data must contain: eid, test_name, test_score, age_at_test, Sex,
#   Edu, instance
#
# USAGE -----------------------------------------------------------------------
#
# Step 1: Source this script to load functions, packages, and data:
#   source("ukb_cog_CV.R")
#
# Step 2: Choose ONE of three execution modes:
#
#   (a) SINGLE TEST (recommended for first run / debugging):
#       result <- run_single_test("BL")
#
#   (b) ALL TESTS, SEQUENTIAL (slow but transparent console output):
#       results <- run_all_sequential()
#
#   (c) ALL TESTS, PARALLEL (faster but suppresses console output):
#       results <- run_all_parallel(workers = 10)
#
# Step 3: Inspect summary statistics:
#   The runner functions print per-test coverage and z-score summaries.
#
# OUTPUT
#   cv_results/
#     |-- cv_results_<test>.csv            # per-test z-scores & percentiles
#     |-- cv_results_all_tests.csv         # combined across all tests
#
# =============================================================================



# 1. SETUP: packages, configuration, and data

suppressPackageStartupMessages({
  library(gamlss)
  library(gamlss.dist)
  library(gamlss.add)      # required for random()
  library(dplyr)
  library(caret)
  library(future.apply)    # parallel backend
  library(parallel)
  library(progressr)       # progress bar for parallel runs
})

options(warn = 1)

# --- Configuration -----------------------------------------------------------
DATA_PATH  <- "path/long_data.RData"
OUTPUT_DIR <- "cv_results"
N_FOLDS    <- 10
CV_SEED    <- 123

# Optimal distribution per test (from prior dist_test results)
# UPDATE THESE after re-running dist_test with random effects
OPTIMAL_DISTRIBUTIONS <- list(
  "BL"   = "ST1",
  "FI"   = "ST2",
  "MP"   = "ST2",
  "NM"   = "BCT",
  "PAL"  = "NO",
  "PMat" = "ST2",
  "PV"   = "BEINF",
  "RT"   = "BCT",
  "SD"   = "ST2",
  "TMA"  = "BCCG",
  "TMB"  = "BCT",
  "TR"   = "NO"
)

# --- Load data ---------------------------------------------------------------
load(DATA_PATH)
cat("Loaded long_data:", nrow(long_data), "rows,",
    length(unique(long_data$eid)), "unique participants\n\n")

# Flag baseline (earliest) observation per participant per test.
# Training data uses baseline only; test data uses ALL instances.
long_data <- long_data %>%
  group_by(eid, test_name) %>%
  arrange(age_at_test) %>%
  mutate(is_baseline = row_number() == 1) %>%
  ungroup()

cat("Baseline observations:", sum(long_data$is_baseline),
    "/", nrow(long_data), "rows\n")

# --- Ensure output directory exists ------------------------------------------
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# 2. HELPER FUNCTIONS

# Map distribution name to family generator + CDF
get_distribution_info <- function(dist_name) {
  dist_map <- list(
    "NO"    = list(family = NO,    pfunc = pNO),
    "BCT"   = list(family = BCT,   pfunc = pBCT),
    "BCCG"  = list(family = BCCG,  pfunc = pBCCG),
    "ST1"   = list(family = ST1,   pfunc = pST1),
    "ST2"   = list(family = ST2,   pfunc = pST2),
    "BE"    = list(family = BE,    pfunc = pBE),
    "BEINF" = list(family = BEINF, pfunc = pBEINF),
    "TF"    = list(family = TF,    pfunc = pTF),
    "SN1"   = list(family = SN1,   pfunc = pSN1)
  )
  if (!dist_name %in% names(dist_map)) stop("Unknown distribution: ", dist_name)
  dist_map[[dist_name]]
}


# 3. CORE ANALYSIS: 10-fold CV for a single cognitive test

cv_single_cognitive <- function(test_data, test_name, dist_name,
                                 output_dir = OUTPUT_DIR, n_folds = N_FOLDS) {

  cat("Testing:", test_name, "(n =", nrow(test_data), ", dist =", dist_name, ")\n")

  dist_info <- get_distribution_info(dist_name)

  # -- Prepare modelling data ------------------------------------------------
  data <- test_data %>%
    dplyr::filter(!is.na(test_score), !is.na(age_at_test),
                  !is.na(Sex), !is.na(Edu)) %>%
    dplyr::transmute(
      eid         = eid,
      instance    = instance,
      is_baseline = is_baseline,
      age         = as.numeric(age_at_test),
      Sex         = factor(Sex),
      Edu         = as.numeric(Edu),
      phenotype   = as.numeric(test_score),
      age_at_test = as.numeric(age_at_test)
    )

  if (nrow(data) == 0) {
    cat("  No valid rows after filtering\n\n")
    return(NULL)
  }
  data$orig_idx <- seq_len(nrow(data))
  data <- as.data.frame(data)

  cat("  Data:", nrow(data), "observations,",
      length(unique(data$eid)), "participants\n")

  # -- Set up 10-fold CV by participant --------------------------------------
  unique_participants  <- unique(data$eid)
  participant_indices  <- seq_along(unique_participants)
  set.seed(CV_SEED)
  folds_idx <- caret::createFolds(participant_indices, k = n_folds,
                                  returnTrain = FALSE)

  # -- Storage ---------------------------------------------------------------
  zscores     <- rep(NA_real_, nrow(data))
  percentiles <- rep(NA_real_, nrow(data))
  conv_status <- integer(n_folds)

  # -- Loop over folds -------------------------------------------------------
  for (k in seq_len(n_folds)) {
    cat("    Fold", k, "/", n_folds)

    test_indices  <- folds_idx[[k]]
    train_indices <- setdiff(participant_indices, test_indices)
    
    test_ids  <- unique_participants[test_indices]
    train_ids <- unique_participants[train_indices]

    data_test  <- droplevels(data[data$eid %in% test_ids, , drop = FALSE])
    data_train <- droplevels(data[data$eid %in% train_ids & data$is_baseline, , drop = FALSE])
    
    cat(" (train:", nrow(data_train), "test:", nrow(data_test), ")")
    if (nrow(data_train) == 0 || nrow(data_test) == 0) {
      cat(" - empty set\n"); next
    }

    # --- Fit GAMLSS model with random intercept for participant ----------
    fit_ok <- FALSE
    mdl    <- NULL

    tryCatch({
      form_mu    <- as.formula("phenotype ~ pb(age) + Sex + pb(Edu)")
      form_sigma <- as.formula("~ pb(age) + Sex + pb(Edu)")

      suppressWarnings({
        mdl <- gamlss::gamlss(
          formula       = form_mu,
          sigma.formula = form_sigma,
          nu.formula    = ~ 1,
          tau.formula   = ~ 1,
          family        = dist_info$family,
          data          = data_train,
          control       = gamlss.control(n.cyc = 400, trace = FALSE, c.crit = 0.01),
          gamlss.env    = environment()
        )
      })

      fit_ok <- isTRUE(mdl$converged) && !any(is.na(fitted(mdl)))

      # Inline data + family into model call so predictAll doesn't need symbol lookup
      if (fit_ok) {
        mdl$call$data       <- data_train
        mdl$call$family     <- dist_info$family
        mdl$call$gamlss.env <- environment()
      }
    }, error = function(e) {
      cat(" ERROR:", substr(e$message, 1, 50))
    })

    if (!fit_ok || is.null(mdl)) { cat(" X (not converged)\n"); next }
    conv_status[k] <- 1

    # --- Predict distribution parameters on test set ---------------------
    params <- tryCatch({
      as.data.frame(predictAll(mdl, newdata = data_test, type = "response"))
    }, error = function(e) {
      cat(" PREDICT ERROR:", substr(e$message, 1, 50)); NULL
    })

    if (is.null(params)) { cat(" X (predict failed)\n"); next }

    # --- Compute percentiles from CDF, convert to z ----------------------
    args <- list(q = as.numeric(data_test$phenotype))
    for (nm in intersect(colnames(params), c("mu", "sigma", "nu", "tau"))) {
      args[[nm]] <- params[[nm]]
    }
    pcts <- do.call(dist_info$pfunc, args)

    eps  <- 1e-10
    pcts <- pmin(pmax(pcts, eps), 1 - eps)
    z    <- qnorm(pcts)
    z[!is.finite(z)]       <- NA_real_
    pcts[!is.finite(pcts)] <- NA_real_

    zscores[data_test$orig_idx]     <- z
    percentiles[data_test$orig_idx] <- 100 * pcts

    cat(" OK\n")
  }

  # -- Merge results back to original test data -----------------------------
  result_data <- test_data %>% dplyr::filter(!is.na(test_score))

  lookup_df <- data.frame(
    eid         = data$eid,
    instance    = data$instance,
    age_at_test = data$age_at_test,
    z_score     = zscores,
    percentile  = percentiles,
    stringsAsFactors = FALSE
  )

  result_data <- dplyr::left_join(result_data, lookup_df,
                                  by = c("eid", "instance", "age_at_test"))

  # -- Save and report ------------------------------------------------------
  filename <- file.path(output_dir, paste0("cv_results_", test_name, ".csv"))
  write.csv(result_data, filename, row.names = FALSE)

  valid_z  <- sum(!is.na(zscores))
  coverage <- round(100 * valid_z / nrow(data), 1)
  cat("  Convergence:", sum(conv_status), "/", n_folds, "folds\n")
  cat("  Coverage:", coverage, "% (", valid_z, "/", nrow(data), "observations)\n")
  cat("  Saved:", filename, "\n\n")

  result_data
}


# 4. WRAPPER: filter to one test, then call core analysis

run_single_test <- function(test_key) {
  cat("=== GAMLSS Cross-Validation:", test_key, "===\n")

  if (!test_key %in% names(OPTIMAL_DISTRIBUTIONS)) {
    stop("Test '", test_key, "' has no optimal distribution defined. ",
         "Available: ", paste(names(OPTIMAL_DISTRIBUTIONS), collapse = ", "))
  }

  dist_name <- OPTIMAL_DISTRIBUTIONS[[test_key]]
  test_data <- long_data %>% dplyr::filter(.data$test_name == test_key)

  if (nrow(test_data) == 0) {
    cat("No data for test:", test_key, "\n")
    return(NULL)
  }

  result <- cv_single_cognitive(test_data, test_key, dist_name)
  cat("=== Complete:", test_key, "===\n\n")
  result
}


# 5. RUNNERS: sequential or parallel execution across all tests

# Run all tests one at a time (best for debugging - full console output)
run_all_sequential <- function() {
  available_tests <- intersect(unique(long_data$test_name),
                               names(OPTIMAL_DISTRIBUTIONS))
  cat("Running", length(available_tests), "tests sequentially...\n\n")

  all_results <- lapply(available_tests, function(tk) run_single_test(tk))
  names(all_results) <- available_tests

  combine_and_summarise(all_results)
}

# Run all tests in parallel (faster, but console output suppressed in workers - silent bugs lay here)
run_all_parallel <- function(workers = 10) {
  available_tests <- intersect(unique(long_data$test_name),
                               names(OPTIMAL_DISTRIBUTIONS))
  cat("Running", length(available_tests), "tests in parallel across",
      workers, "workers...\n\n")

  plan(multisession, workers = workers)
  progressr::handlers(global = TRUE)
  progressr::handlers("txtprogressbar")

  all_results <- progressr::with_progress({
    p <- progressr::progressor(along = available_tests)
    future_lapply(available_tests, function(tk) {
      # Re-attach packages inside each worker session
      library(gamlss); library(gamlss.dist); library(gamlss.add)
      library(dplyr);  library(caret)
      res <- run_single_test(tk)
      p(sprintf("done: %s", tk))
      res
    }, future.seed = TRUE)
  })
  names(all_results) <- available_tests

  combine_and_summarise(all_results)
}


# 6. UTILITIES
                        
# Combine per-test results, save combined file, print summary
combine_and_summarise <- function(all_results) {
  all_results <- all_results[!sapply(all_results, is.null)]
  if (length(all_results) == 0) {
    cat("No results to combine.\n")
    return(invisible(NULL))
  }

  combined_data <- bind_rows(all_results)
  combined_file <- file.path(OUTPUT_DIR, "cv_results_all_tests.csv")
  write.csv(combined_data, combined_file, row.names = FALSE)
  cat("\nAll results saved to:", combined_file, "\n")

  summary_stats <- combined_data %>%
    group_by(test_name) %>%
    summarise(
      n_total      = n(),
      n_with_z     = sum(!is.na(z_score)),
      coverage_pct = round(100 * n_with_z / n_total, 1),
      mean_z       = round(mean(z_score, na.rm = TRUE), 3),
      sd_z         = round(sd(z_score, na.rm = TRUE), 3),
      .groups      = "drop"
    )
  cat("\n=== Summary Statistics ===\n")
  print(summary_stats)

  combined_data
}

show_available_tests <- function() {
  available_tests <- intersect(unique(long_data$test_name),
                               names(OPTIMAL_DISTRIBUTIONS))
  cat("Available cognitive tests:\n")
  for (i in seq_along(available_tests)) {
    tk <- available_tests[i]
    n  <- nrow(long_data %>% dplyr::filter(.data$test_name == tk))
    cat(sprintf("%2d. %-6s [%s] (n = %5d)\n",
                i, tk, OPTIMAL_DISTRIBUTIONS[[tk]], n))
  }
  invisible(available_tests)
}


# -----------------------------------------------------------------------------
# Script loaded -- see USAGE block at top of file for execution commands
# -----------------------------------------------------------------------------
cat("Script loaded\n")
cat("Quick start:\n")
cat("Run single test: result <- run_single_test('BL')\n")
cat("Run all in series: results <- run_all_sequential()\n")
cat("Run all in parallel: results <- run_all_parallel(workers = 10)\n")
cat("List available tests: show_available_tests()\n")
