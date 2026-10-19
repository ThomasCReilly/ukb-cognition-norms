
# GAMLSS Distribution Testing for UK Biobank Cognitive Tests
#
# This is the script for testing the normative distribution families to call in GAMLSS. 
#For each cognitive test, it fits GAMLSS normative models with 9 candidate
# distribution families, then ranks them by AIC, saves diagnostic plots
# and ranked model comparison tables to disk.
#
# MODEL STRUCTURE
#   mu     ~ pb(age_at_test) + Sex + pb(Edu) + random(as.factor(eid))
#   sigma  ~ pb(age_at_test) + Sex + pb(Edu)
#   nu     ~ 1
#   tau    ~ 1
#
#   The random() term is a participant-level random intercept that accounts for
#   within-participant correlation across repeated cognitive assessments.
#
# DATA REQUIREMENTS
#   long_data must contain: eid, test_name, test_score, age_at_test, Sex, Edu
#   (loaded automatically from DATA_PATH below)
#
# USAGE -----------------------------------------------------------------------
#
# Step 1: Source this script to load functions, packages, and data:
#   source("dist_test.R")
#
# Step 2: Choose ONE of three execution modes:
#
#   (a) SINGLE TEST (recommended for first run / debugging):
#       result <- run_distribution_test(long_data, "BL")
#
#   (b) ALL TESTS, SEQUENTIAL (slow but transparent console output):
#       results_list <- run_all_sequential(long_data)
#
#   (c) ALL TESTS, PARALLEL (faster but suppresses console output):
#       results_list <- run_all_parallel(long_data, workers = 10)
#
# Step 3: Inspect best distribution per test:
#   best_by_test <- sapply(results_list, function(x) x$best_distribution)
#
# OUTPUT
#   gamlss_results/<test_name>/
#     |-- <test>_model_comparison.csv     # all distributions tested
#     |-- <test>_ranked_models.csv        # successful models, ranked by AIC
#     |-- <test>_residual_stats.csv       # residual diagnostics
#     |-- <test>_results.RData            # full fitted models
#     |-- residual_plots/<dist>_residuals.png  # 2x2 diagnostic panel per dist
#


# 1. SETUP: packages, configuration, and data

suppressPackageStartupMessages({
  library(gamlss)
  library(gamlss.dist)
  library(gamlss.add)      # required for re()
  library(dplyr)
  library(moments)
  library(splines)
  library(future.apply)    # parallel backend
  library(parallel)
})

options(warn = 1)

# --- Configuration -----------------------------------------------------------
DATA_PATH     <- "path/long_data.RData" #replace with data path
OUTPUT_DIR    <- "gamlss_results"
DISTRIBUTIONS <- c("ST1", "ST2", "BCCG", "TF", "SN1", "NO", "BE", "BEINF", "BCT")
MIN_N         <- 50

# --- Load data ---------------------------------------------------------------
load(DATA_PATH)
cat("Loaded long_data:", nrow(long_data), "rows,",
    length(unique(long_data$eid)), "unique participants\n\n")



# 2. HELPER FUNCTIONS

# Create directory for a single test's output
create_test_directory <- function(base_dir, test_name) {
  test_dir <- file.path(base_dir, test_name)
  dir.create(test_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(test_dir, "residual_plots"), showWarnings = FALSE)
  test_dir
}

# Smooth mean +/- SD over a continuous covariate (used for age diagnostic plot)
draw_mean_sd_smooth <- function(x, r, xlab, main, xlim_low = NULL, span = 0.30) {
  ok <- is.finite(x) & is.finite(r)
  x <- x[ok]; r <- r[ok]
  if (length(unique(x)) < 10L || length(x) < 30L) {
    plot.new(); title(main = paste(main, "- Non-calculable")); return(invisible())
  }
  o <- order(x); x <- x[o]; r <- r[o]
  
  m_fit <- try(loess(r ~ x, span = span, degree = 1, family = "symmetric"), silent = TRUE)
  v_fit <- try(loess(I(r^2) ~ x, span = span, degree = 1, family = "symmetric"), silent = TRUE)
  if (inherits(m_fit, "try-error") || inherits(v_fit, "try-error")) {
    plot.new(); title(main = paste(main, "- Non-calculable")); return(invisible())
  }
  
  xr <- range(x, na.rm = TRUE)
  if (!is.null(xlim_low)) xr[1] <- min(xlim_low, xr[1])
  xg <- seq(xr[1], xr[2], length.out = 200)
  
  mu  <- as.numeric(predict(m_fit, newdata = data.frame(x = xg)))
  m2  <- as.numeric(predict(v_fit, newdata = data.frame(x = xg)))
  sdv <- sqrt(pmax(m2 - mu^2, 0))
  
  ylo <- mu - sdv; yhi <- mu + sdv
  yr  <- range(c(ylo, yhi, 0), na.rm = TRUE)
  pad <- 0.06 * diff(yr)
  ylim <- c(yr[1] - pad, yr[2] + pad)
  
  old <- par(xaxs = "r"); on.exit(par(old), add = TRUE)
  plot(xg, mu, type = "n",
       xlab = xlab, ylab = "Residual mean +/- SD",
       main = main, xlim = xr, ylim = ylim)
  polygon(c(xg, rev(xg)), c(ylo, rev(yhi)),
          col = rgb(0, 0, 0, 0.12), border = NA)
  lines(xg, mu, lwd = 2)
  abline(h = 0, col = "red", lty = 2)
}


# 3. CORE ANALYSIS: fit & evaluate all distributions for a single test


test_single_cognitive <- function(test_data, test_name, output_dir = OUTPUT_DIR) {
  
  cat("Testing:", test_name, "(n =", nrow(test_data), ")\n")
  
  test_dir <- create_test_directory(output_dir, test_name)
  
  # Storage
  results <- data.frame(
    distribution = character(),
    aic          = numeric(),
    bic          = numeric(),
    converged    = logical(),
    error_msg    = character(),
    stringsAsFactors = FALSE
  )
  residual_stats <- data.frame(
    distribution = character(),
    mean         = numeric(),
    sd           = numeric(),
    skewness     = numeric(),
    kurtosis     = numeric(),
    n_infinite   = integer(),
    note         = character(),
    stringsAsFactors = FALSE
  )
  models <- list()
  
  age_vec <- test_data$age_at_test
  edu_vec <- test_data$Edu
  
  # -- Loop over candidate distributions --------------------------------------
  for (dist in DISTRIBUTIONS) {
    cat("  ", dist, "... ")
    
    result_row <- list(
      distribution = dist,
      aic          = NA,
      bic          = NA,
      converged    = FALSE,
      error_msg    = ""
    )
    
    tryCatch({
      
      # --- Fit GAMLSS model with random intercept for participant ----------
      # Build formulas with explicit environment binding

      form_mu <- as.formula("test_score ~ pb(age_at_test) + Sex + pb(Edu) + random(as.factor(eid))")
      form_sigma <- as.formula("~ pb(age_at_test) + Sex + pb(Edu)")
      
      suppressWarnings({
        model <- gamlss(
          formula       = form_mu,
          sigma.formula = form_sigma,
          nu.formula    = ~ 1,
          tau.formula   = ~ 1,
          family        = get(dist)(),
          data          = test_data,
          control       = gamlss.control(trace = FALSE, n.cyc = 200)
        )
      })
      
      if (isTRUE(model$converged) && !any(is.na(fitted(model)))) {
        result_row$aic       <- AIC(model)
        result_row$bic       <- BIC(model)
        result_row$converged <- TRUE
        models[[dist]]       <- model
        
        # --- Residuals (compute once, reuse) -----------------------------
        qres_try <- try(residuals(model, what = "z"), silent = TRUE)
        qres_raw <- if (inherits(qres_try, "try-error")) {
          rep(NA_real_, nrow(test_data))
        } else qres_try
        qrf   <- qres_raw[is.finite(qres_raw)]
        n_inf <- sum(!is.finite(qres_raw))
        
        # --- Residual statistics -----------------------------------------
        if (length(qrf) >= 3) {
          residual_stats <- rbind(residual_stats, data.frame(
            distribution = dist,
            mean         = round(mean(qrf), 4),
            sd           = round(sd(qrf), 4),
            skewness     = round(skewness(qrf), 4),
            kurtosis     = round(kurtosis(qrf), 4),  # excess kurtosis
            n_infinite   = n_inf,
            note         = "ok",
            stringsAsFactors = FALSE
          ))
        } else {
          residual_stats <- rbind(residual_stats, data.frame(
            distribution = dist,
            mean = NA_real_, sd = NA_real_,
            skewness = NA_real_, kurtosis = NA_real_,
            n_infinite = n_inf,
            note = "Residuals not calculable",
            stringsAsFactors = FALSE
          ))
        }
        
        # --- 2x2 diagnostic plots (QQ, hist, age, edu) -------------------
        res_file <- file.path(test_dir, "residual_plots",
                              paste0(dist, "_residuals.png"))
        blank <- function(t) { plot.new(); title(main = t) }
        
        png(res_file, width = 15, height = 12, units = "in", res = 600, bg = "white")
        par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))
        
        ## (1) QQ plot
        if (length(qrf) >= 3) {
          qqnorm(qrf, main = paste("QQ Plot:", dist),
                 xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
                 col = "black", pch = 20, cex = 0.5)
          qqline(qrf, col = "red", lty = 2)
        } else blank(paste("QQ Plot:", dist, "- Non-calculable"))
        
        ## (2) Histogram + density + N(0,1) overlay
        if (length(qrf) >= 3) {
          h <- hist(qrf, breaks = 30, plot = FALSE)
          d <- density(qrf, na.rm = TRUE)
          ylim_hist <- c(0, max(h$density, d$y, na.rm = TRUE) * 1.025)
          plot(h, freq = FALSE,
               main = paste("Residual Histogram:", dist),
               xlab = "Standardized residuals",
               col = "grey85", border = "black", ylim = ylim_hist)
          lines(d, col = "black", lwd = 2)
          curve(dnorm(x, 0, 1), add = TRUE, col = "red", lty = 2)
        } else blank(paste("Residual Histogram:", dist, "- Non-calculable"))
        
        ## (3) Mean +/- SD vs age (LOESS)
        if (length(qres_raw) == length(age_vec)) {
          draw_mean_sd_smooth(age_vec, qres_raw,
                              xlab = "Age",
                              main = paste(dist, "- Mean +/- SD by Age"),
                              xlim_low = 35, span = 0.30)
        } else blank(paste(dist, "- Mean +/- SD by Age - Non-calculable"))
        
        ## (4) Mean +/- SD vs education (discrete)
        ok <- is.finite(edu_vec) & is.finite(qres_raw)
        edu <- edu_vec[ok]; r <- qres_raw[ok]
        if (length(unique(edu)) >= 2 && length(r) >= 3) {
          edu_levels <- sort(unique(edu))
          mns <- tapply(r, edu, mean, na.rm = TRUE)
          sds <- tapply(r, edu, sd,   na.rm = TRUE)
          sds[is.na(sds)] <- 0  # singleton groups
          ylo <- mns - sds; yhi <- mns + sds
          ylim <- range(c(mns, ylo, yhi, 0), na.rm = TRUE)
          plot(edu_levels, mns, type = "b", pch = 19, lwd = 2,
               xlab = "Education Rank", ylab = "Residual mean +/- SD",
               main = paste(dist, "- Mean +/- SD by Edu"), ylim = ylim)
          arrows(edu_levels, ylo, edu_levels, yhi,
                 angle = 90, code = 3, length = 0.05, col = "grey50")
          abline(h = 0, col = "red", lty = 2)
        } else blank(paste(dist, "- Mean +/- SD by Edu - Non-calculable"))
        
        par(mfrow = c(1, 1))
        dev.off()
        
        cat("AIC:", round(result_row$aic, 1), "\n")
        
      } else {
        result_row$error_msg <- "No convergence"
        cat("FAILED\n")
      }
      
    }, error = function(e) {
      result_row$error_msg <<- substr(e$message, 1, 50)
      cat("ERROR:", e$message, "\n")
    })
    
    results <- rbind(results, data.frame(result_row, stringsAsFactors = FALSE))
  }
  
  # -- Rank successful models by AIC ------------------------------------------
  successful <- results[results$converged, ]
  if (nrow(successful) > 0) {
    successful <- successful[order(successful$aic), ]
    successful$delta_aic <- successful$aic - min(successful$aic)
    successful$rank      <- seq_len(nrow(successful))
    successful$support   <- ifelse(successful$delta_aic < 2, "Strong",
                                   ifelse(successful$delta_aic < 7, "Some", "Little"))
  }
  
  # -- Save outputs -----------------------------------------------------------
  write.csv(results,        file.path(test_dir, paste0(test_name, "_model_comparison.csv")), row.names = FALSE)
  write.csv(successful,     file.path(test_dir, paste0(test_name, "_ranked_models.csv")),    row.names = FALSE)
  write.csv(residual_stats, file.path(test_dir, paste0(test_name, "_residual_stats.csv")),   row.names = FALSE)
  save(models, results, residual_stats,
       file = file.path(test_dir, paste0(test_name, "_results.RData")))
  
  cat("Results saved to:", test_dir, "\n\n")
  
  list(
    results           = results,
    successful        = successful,
    residual_stats    = residual_stats,
    models            = models,
    best_distribution = if (nrow(successful) > 0) successful$distribution[1] else NA,
    directory         = test_dir
  )
}


# 4. WRAPPER: filter to one test, then call core analysis


run_distribution_test <- function(long_data, test_key, min_n = MIN_N) {
  cat("=== GAMLSS Distribution Testing:", test_key, "===\n")
  
  test_data <- long_data %>% dplyr::filter(.data$test_name == test_key)
  if (nrow(test_data) < min_n) {
    cat("Insufficient data (n =", nrow(test_data), "). Need at least", min_n, "\n")
    return(NULL)
  }
  
  result <- test_single_cognitive(test_data, test_key)
  cat("=== Complete:", test_key, "===\n\n")
  result
}



# 5. RUNNERS: sequential or parallel execution across all tests

# Run all tests one at a time (best for debugging - full console output)
run_all_sequential <- function(long_data) {
  tests <- unique(long_data$test_name)
  cat("Running", length(tests), "tests sequentially...\n\n")
  
  results_list <- lapply(tests, function(tk) run_distribution_test(long_data, tk))
  names(results_list) <- tests
  results_list
}

# Run all tests in parallel (faster, but console output suppressed in workers)
run_all_parallel <- function(long_data, workers = 10) {
  tests <- unique(long_data$test_name)
  cat("Running", length(tests), "tests in parallel across", workers, "workers...\n\n")
  
  plan(multisession, workers = workers)
  progressr::handlers(global = TRUE)
  progressr::handlers("txtprogressbar")
  
  results_list <- progressr::with_progress({
    p <- progressr::progressor(along = tests)
    future_lapply( tests, function(tk) {
      # Re-attach packages inside each worker session
      library(gamlss); library(gamlss.dist); library(gamlss.add)
      library(dplyr);  library(moments);     library(splines)
      run_distribution_test(long_data, tk)
      res <- run_distribution_test(long_data, tk)
      p(sprintf("done: %s", tk))
      res
    }, future.seed = TRUE)
  })
  names(results_list) <- tests
  results_list
}


# 6. UTILITIES

show_available_tests <- function(long_data) {
  tests <- unique(long_data$test_name)
  cat("Available cognitive tests:\n")
  for (i in seq_along(tests)) {
    n <- nrow(long_data %>% dplyr::filter(.data$test_name == tests[i]))
    cat(sprintf("%2d. %-15s (n = %5d)\n", i, tests[i], n))
  }
  invisible(tests)
}


# -----------------------------------------------------------------------------
# Script loaded -- see USAGE block at top of file for execution commands
# -----------------------------------------------------------------------------
cat("Script loaded\n")
cat("Quick start:\n")
cat("Run single test: result <- run_distribution_test(long_data, 'BL')\n")
cat("Run all in series: results_list <- run_all_sequential(long_data)\n")
cat("Run all in parallel: results_list <- run_all_parallel(long_data, workers = 10)\n")
