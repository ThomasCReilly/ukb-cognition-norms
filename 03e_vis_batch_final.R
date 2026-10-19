#!/usr/bin/env Rscript

# This script is used to generate visualisations of the GAMLSS models by
# showing the predicted test scores accross age and education. 
# Visualise normative model surfaces from fitted GAMLSS models.
# For each cognitive test and sex:
#   1. Render a 3D mesh surface (median coloured by band width, +/- SD shells)
#   2. Plot a cross-sectional slice at a fixed age (score ~ education)
#   3. Plot a cross-sectional slice at a fixed education level (score ~ age)
#
# Bands are CENTILE-based: "+/- k SD" means the pnorm(+/-k) centiles of the
# fitted distribution (all of mu, sigma, nu, tau), i.e. z = +/-k. This is
# correct for any GAMLSS family, unlike mu +/- k*sigma, which only equals
# mean +/- k SD for the normal (NO) family.

#To run for all - source("vis_batch_final.R")
#To run for one test - run_for_test("TMB")

# SETUP

suppressPackageStartupMessages({
  library(gamlss)
  library(dplyr)
  library(ggplot2)
  library(svglite)
  library(plotly)
  library(htmlwidgets)
})

# CONFIGURATION
# EDU_RANGE is descending by design (3D axis orientation).

MODELS_ROOT <- Sys.getenv("MODELS_ROOT", "path/gamlss_results_baseline")
LONG_DATA   <- Sys.getenv("LONG_DATA_PATH", "path/long_data.RData")
OUTDIR      <- Sys.getenv("OUTDIR", "surfaces")

SEXES       <- c("Male", "Female")
AGE_RANGE   <- c(35, 90)
EDU_RANGE   <- c(6, 1)
N_POINTS    <- 200   # grid points per axis (was 700 = 490k cells per surface)
SD_LEVELS   <- 1.5   # z-level(s) for bands; slices use max(SD_LEVELS)

CROSS_SECTION_AGE <- 60
CROSS_SECTION_EDU <- 5

# Observed-data overlay on cross-sections
AGE_WINDOW       <- 2.5    # years either side of CROSS_SECTION_AGE
MAX_POINTS       <- 3000   # subsample cap per panel (keeps SVGs light)
POINT_SEED       <- 2024
RASTERISE_POINTS <- TRUE   # rasterise points inside SVG if ggrastr installed

# Text sizes and 3D canvas
SURF_FONT      <- "Arial"
SURF_TITLE_PT  <- 22   # 3D title
SURF_AXIS_PT   <- 18   # 3D axis titles
SURF_TICK_PT   <- 15   # 3D tick labels + colourbar
SURF_WIDTH     <- 760  # 3D layout size in px; static export uses the same
SURF_HEIGHT    <- 670  #   size and multiplies resolution by `scale`
# Camera. The scene is clipped to its own box (`domain`), so to move the mesh
# sideways within that box, shift eye and centre together: +x/-y moves the
# mesh right on screen. Eye distance from centre controls zoom.
SURF_SHIFT     <- 0.04   # sideways shift (0 = centred)
SURF_EYE       <- list(x = 1.33 + SURF_SHIFT, y = 1.33 - SURF_SHIFT, z = 0.62)
SURF_CENTER    <- list(x = SURF_SHIFT, y = -SURF_SHIFT, z = -0.24)
SLICE_BASE_PT  <- 18   # ggplot base_size for cross-sections

STRICT_FAMILY_LOOKUP <- FALSE

OPTIMAL_DISTRIBUTIONS <- list(
  "BL"   = "ST1",  "FI"   = "ST2",  "MP"   = "ST2",  "NM"   = "BCT",
  "PAL"  = "NO",   "PMat" = "ST2",  "PV"   = "BEINF","RT"   = "BCT",
  "SD"   = "ST2",  "TMA"  = "BCCG", "TMB"  = "BCT",  "TR"   = "NO"
)

SCORE_RANGES <- list(
  "BL"   = c(5, 40),   "FI"   = c(0, 13),  "MP"   = c(0, 15),
  "NM"   = c(4, 8),    "PAL"  = c(0, 10),  "PMat" = c(0, 10),
  "PV"   = c(0, 1),    "RT"   = c(250, 750),"SD"   = c(0, 40),
  "TMA"  = c(10, 40),  "TMB"  = c(20, 125), "TR"   = c(0, 18)
)

SCORE_UNITS <- list(
  "BL"   = "n correct",  "FI"   = "n correct",  "MP"  = "n correct",
  "NM"   = "n correct",  "PAL"  = "n correct",  "PMat"= "n errors",
  "PV"   = "index",      "RT"   = "ms",         "SD"  = "n correct",
  "TMA"  = "s",          "TMB"  = "s",          "TR"  = "n correct"
)

# LOAD LONG DATA
# Must contain: eid, Sex, Edu, test_name, instance, age_at_test, test_score.

ld_env <- new.env(parent = emptyenv())
objs   <- load(LONG_DATA, envir = ld_env)
stopifnot("long_data" %in% objs)
long_data <- ld_env$long_data
stopifnot(all(
  c("eid", "Sex", "Edu", "test_name", "instance", "age_at_test", "test_score")
  %in% names(long_data)
))

# Keep baseline (earliest) observation per participant per test.
filter_to_baseline <- function(data) {
  data %>%
    dplyr::group_by(eid, test_name) %>%
    dplyr::arrange(age_at_test) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
}

long_data <- filter_to_baseline(long_data)

# SMALL HELPERS

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
fmt_n  <- function(n) format(n, big.mark = ",", scientific = FALSE)

score_label <- function(test_name) {
  u <- SCORE_UNITS[[test_name]]
  if (!is.null(u) && nzchar(u)) paste0("Test Score (", u, ")") else "Test Score"
}

# Subset rows to one sex, handling numeric (1 = Male) or character coding.
filter_sex <- function(d, sex) {
  if (is.null(sex) || !"Sex" %in% names(d)) return(d)
  keep <- if (is.numeric(d$Sex)) {
    d$Sex == ifelse(tolower(sex) == "male", 1, 0)
  } else {
    tolower(trimws(as.character(d$Sex))) == tolower(sex)
  }
  d[which(keep), , drop = FALSE]
}

# MODEL I/O HELPERS

models_file_for <- function(test) {
  file.path(MODELS_ROOT, test, sprintf("%s_results.RData", test))
}

load_models_for_test <- function(test) {
  rfile <- models_file_for(test)
  if (!file.exists(rfile)) stop("Missing model file: ", rfile)
  e  <- new.env(parent = emptyenv())
  ob <- load(rfile, envir = e)
  if (!"models" %in% ob) stop("'models' object not in ", rfile)
  e$models
}

# MODEL SELECTION
# Exact list-name match first, then fall back to each model's family slot.

pick_model <- function(models, fam, sex = NULL) {
  nm <- names(models)
  
  if (!is.null(nm) && fam %in% nm) return(models[[fam]])
  
  if (isTRUE(STRICT_FAMILY_LOOKUP)) {
    stop("Family '", fam, "' not found in names(models). Have: ",
         paste0(nm, collapse = ", "))
  }
  
  same_fam <- function(m) {
    ff <- tryCatch(m$family$family, error = function(e) NA_character_)
    !is.na(ff[1]) && toupper(ff[1]) == toupper(fam)
  }
  idx <- which(vapply(models, same_fam, logical(1)))
  
  if (length(idx) == 1) return(models[[idx]])
  if (length(idx) > 1 && !is.null(nm) && !is.null(sex)) {
    pref <- idx[grep(tolower(sex), tolower(nm)[idx])]
    if (length(pref) >= 1) return(models[[pref[1]]])
  }
  if (length(idx) >= 1) return(models[[idx[1]]])
  
  stop("No model of family '", fam, "' found. names(models): ",
       paste0(nm, collapse = ", "))
}

# OBSERVED SCORE RANGE
# Pools across sexes by default so male/female plots share a scale.

score_z_range <- function(df, test, sex = NULL, use_all_sexes = TRUE,
                          pad_abs = 0, min_span = 0) {
  d <- dplyr::filter(df, .data$test_name == !!test)
  if (!use_all_sexes) d <- filter_sex(d, sex)
  
  if (!nrow(d)) return(c(0, 1))
  v <- as_num(d$test_score)
  v <- v[is.finite(v)]
  if (!length(v)) return(c(0, 1))
  
  lo <- min(v); hi <- max(v)
  if (lo == hi) {
    span <- max(1e-6, min_span)
    return(c(lo - span / 2, hi + span / 2))
  }
  if (min_span > 0 && (hi - lo) < min_span) {
    mid <- (lo + hi) / 2
    lo  <- mid - min_span / 2
    hi  <- mid + min_span / 2
  }
  c(lo - pad_abs, hi + pad_abs)
}

# PREPARE MODEL FOR PREDICTION
# Re-attach source data and family so prediction on new grids is reliable.

prepare_model <- function(mdl, fam, test_data) {
  mdl$call$data   <- test_data
  mdl$call$family <- as.gamlss.family(fam)
  mdl
}

sex_value <- function(sex, long_data) {
  if (is.numeric(long_data$Sex)) {
    ifelse(sex == "Male", 1, 0)
  } else {
    factor(sex, levels = c("Female", "Male"))
  }
}

# CENTILE PREDICTION
# predict_params(): all distribution parameters on the response scale.
# centile(): quantile at z (e.g. z = 0 -> median, z = 1.5 -> 93.3rd centile)
# using the family's own q-function, so skew/kurtosis are respected.

predict_params <- function(mdl, newdata, test_data) {
  pnames <- intersect(c("mu", "sigma", "nu", "tau"), mdl$parameters)
  setNames(
    lapply(pnames, function(w)
      as.numeric(predict(mdl, newdata = newdata, what = w, type = "response"))),
    pnames
  )
}

centile <- function(pars, fam, z) {
  qfun <- get(paste0("q", fam), mode = "function")
  n    <- length(pars$mu)
  do.call(qfun, c(list(p = rep(pnorm(z), n)), pars))
}

# 3D SURFACE VISUALISATION
# Median surface over age x education, coloured by the width of the
# +/- k SD band (a skew-robust measure of spread), with translucent shells.

visualise_surface <- function(test_name, long_data, models,
                              age_range  = AGE_RANGE,
                              edu_range  = EDU_RANGE,
                              z_range    = NULL,
                              n_points   = N_POINTS,
                              sex        = "Male",
                              sd_levels  = SD_LEVELS) {
  
  fam <- OPTIMAL_DISTRIBUTIONS[[test_name]]
  stopifnot(!is.null(fam))
  
  test_data <- dplyr::filter(long_data, .data$test_name == !!test_name)
  mdl       <- prepare_model(pick_model(models, fam, sex), fam, test_data)
  
  edu_seq <- seq(edu_range[1], edu_range[2], length.out = n_points)
  age_seq <- seq(age_range[1], age_range[2], length.out = n_points)
  
  pred_grid <- expand.grid(
    Edu         = edu_seq,
    age_at_test = age_seq,
    Sex         = sex_value(sex, long_data)
  )
  pars <- predict_params(mdl, pred_grid, test_data)
  
  # expand.grid varies Edu fastest, so each run of nx values is one age.
  # Fill row-wise: rows = age (plotly y), columns = edu (plotly x).
  # Correct for any nx/ny, not just nx == ny.
  nx <- length(edu_seq); ny <- length(age_seq)
  to_mat <- function(v) matrix(v, nrow = ny, ncol = nx, byrow = TRUE)
  
  med_mat <- to_mat(centile(pars, fam, 0))
  bands <- lapply(sd_levels, function(k) list(
    k  = k,
    hi = to_mat(centile(pars, fam,  k)),
    lo = to_mat(centile(pars, fam, -k))
  ))
  k_col     <- max(sd_levels)
  b_col     <- bands[[which.max(sd_levels)]]
  width_mat <- b_col$hi - b_col$lo
  
  if (is.null(z_range)) {
    z_range <- score_z_range(long_data, test_name, use_all_sexes = TRUE)
  }
  
  ax_title <- function(txt) {
    list(text = paste0("<b>", txt, "</b>"), font = list(size = SURF_AXIS_PT))
  }
  tick <- list(size = SURF_TICK_PT)
  
  # Mostly ambient light: keeps viridis colours true instead of shading
  # parts of the surface dark (which reads as a different band width).
  flat_light <- list(ambient = 0.85, diffuse = 0.25, specular = 0.05,
                     roughness = 0.9, fresnel = 0.1)
  
  p <- plot_ly() |>
    add_surface(
      x = ~edu_seq, y = ~age_seq, z = ~med_mat,
      surfacecolor = width_mat, colorscale = "Viridis",
      lighting = flat_light,
      colorbar = list(
        title = list(text = "<b>Spread</b>",
                     font = list(size = SURF_TICK_PT), side = "top"),
        tickfont = tick, len = 0.6, thickness = 16,
        x = 0.88, xanchor = "left", y = 0.5, yanchor = "middle"
      ),
      name = "Median", showlegend = FALSE
    )
  
  shell_col <- "#9E9E9E"
  for (b in bands) {
    upper <- pmin(b$hi, z_range[2])
    lower <- pmax(b$lo, z_range[1])
    p <- p |>
      add_surface(x = ~edu_seq, y = ~age_seq, z = upper,
                  colorscale = list(c(0, shell_col), c(1, shell_col)),
                  opacity = 0.22, showscale = FALSE, lighting = flat_light,
                  showlegend = FALSE, name = sprintf("+%s SD", b$k)) |>
      add_surface(x = ~edu_seq, y = ~age_seq, z = lower,
                  colorscale = list(c(0, shell_col), c(1, shell_col)),
                  opacity = 0.22, showscale = FALSE, lighting = flat_light,
                  showlegend = FALSE, name = sprintf("\u2212%s SD", b$k))
  }
  
  p |> layout(
    font  = list(family = SURF_FONT, color = "#333333"),
    title = list(
      text = sprintf("<b>%s \u2014 %s</b> (%s)", test_name, sex, fam),
      font = list(size = SURF_TITLE_PT),
      x = 0.5, xanchor = "center", y = 0.97, yanchor = "top"
    ),
    scene = list(
      domain = list(x = c(0, 0.88), y = c(0, 0.95)),
      xaxis = list(title = ax_title("Education (rank)"), tickfont = tick,
                   range = edu_range, autorange = FALSE, dtick = 1),
      yaxis = list(title = ax_title("Age (years)"), tickfont = tick,
                   range = age_range, autorange = FALSE, dtick = 10),
      zaxis = list(title = ax_title(score_label(test_name)), tickfont = tick,
                   range = z_range, autorange = FALSE, nticks = 6),
      aspectmode  = "manual",
      aspectratio = list(x = 1, y = 1, z = 1),
      camera = list(
        eye    = SURF_EYE,
        center = SURF_CENTER,
        up     = list(x = 0, y = 0, z = 1)
      )
    ),
    showlegend = FALSE,
    width  = SURF_WIDTH, height = SURF_HEIGHT,
    margin = list(l = 0, r = 0, t = 40, b = 0)
  )
}

# CROSS-SECTIONAL PLOTS
# 2D slices through the normative surface (ggplot2 -> SVG):
#   - slice_by_age: fix age, score ~ education
#   - slice_by_edu: fix education, score ~ age
# Each shows observed baseline data (in the panel colour), the centile band, the median,
# and direct labels for the band edges on a right-hand axis.

theme_pub <- function(base_size = SLICE_BASE_PT) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
      panel.border     = element_rect(colour = "grey30", fill = NA, linewidth = 0.4),
      axis.ticks       = element_line(colour = "grey30", linewidth = 0.4),
      axis.line        = element_blank(),
      axis.title       = element_text(face = "bold"),
      axis.text        = element_text(colour = "grey20"),
      plot.title       = element_text(face = "bold", hjust = 0.5, size = rel(1.15)),
      plot.subtitle    = element_text(hjust = 0.5, colour = "grey40"),
      plot.caption     = element_text(hjust = 0, colour = "grey40", size = rel(0.6),
                                      lineheight = 1.1),
      legend.position  = "none"
    )
}

# Observed baseline scores for the slice, same sex, subsampled if large.
# window > 0 : keep |fixed_var - fixed_val| <= window (continuous, e.g. age)
# window == 0: keep fixed_var == fixed_val (discrete, e.g. education rank)
observed_points <- function(test_data, sex, x_var, fixed_var, fixed_val,
                            window = 0, max_n = MAX_POINTS, seed = POINT_SEED) {
  d  <- filter_sex(test_data, sex)
  fx <- as_num(d[[fixed_var]])
  keep <- if (window > 0) abs(fx - fixed_val) <= window else abs(fx - fixed_val) < 1e-8
  d <- d[which(keep), , drop = FALSE]
  
  out <- data.frame(x = as_num(d[[x_var]]), y = as_num(d$test_score))
  out <- out[is.finite(out$x) & is.finite(out$y), , drop = FALSE]
  
  n_total <- nrow(out)
  if (n_total > max_n) {
    set.seed(seed)
    out <- out[sample.int(n_total, max_n), , drop = FALSE]
  }
  attr(out, "n_total") <- n_total
  out
}

points_layer <- function(pts, jitter_w, jitter_h, col) {
  if (!nrow(pts)) return(NULL)
  lyr <- geom_point(
    data = pts, aes(x = x, y = y), inherit.aes = FALSE,
    position = position_jitter(width = jitter_w, height = jitter_h, seed = POINT_SEED),
    shape = 16, size = 0.8, colour = col, alpha = 0.2
  )
  if (isTRUE(RASTERISE_POINTS) && requireNamespace("ggrastr", quietly = TRUE)) {
    lyr <- ggrastr::rasterise(lyr, dpi = 300)
  }
  lyr
}

# Label positions at the right-hand end of the curves, nudged apart if close.
band_label_positions <- function(df, k, y_range) {
  i    <- which.max(df$x)
  y    <- c(df$hi[i], df$med[i], df$lo[i])
  span <- diff(y_range)
  inset <- 0.02 * span
  y   <- pmin(pmax(y, y_range[1] + inset), y_range[2] - inset)
  gap <- 0.08 * span
  y[1] <- max(y[1], y[2] + gap)
  y[3] <- min(y[3], y[2] - gap)
  y   <- pmin(pmax(y, y_range[1] + inset), y_range[2] - inset)
  data.frame(
    y   = y,
    lab = c(sprintf("+%s SD", k), "Median", sprintf("\u2212%s SD", k))
  )
}

plot_slice <- function(test_name, long_data, models, sex,
                       x_var, x_seq, fixed_var, fixed_val, window,
                       col, x_lab, xlim, subtitle, point_desc,
                       sd_levels, y_range, jitter_w) {
  
  fam       <- OPTIMAL_DISTRIBUTIONS[[test_name]]
  test_data <- dplyr::filter(long_data, .data$test_name == !!test_name)
  mdl       <- prepare_model(pick_model(models, fam, sex), fam, test_data)
  k         <- max(sd_levels)
  
  grid <- data.frame(
    Edu         = if (x_var == "Edu")         x_seq else fixed_val,
    age_at_test = if (x_var == "age_at_test") x_seq else fixed_val,
    Sex         = sex_value(sex, long_data)
  )
  pars <- predict_params(mdl, grid, test_data)
  df <- data.frame(
    x   = x_seq,
    med = centile(pars, fam, 0),
    lo  = centile(pars, fam, -k),
    hi  = centile(pars, fam,  k)
  )
  
  if (is.null(y_range)) {
    y_range <- score_z_range(long_data, test_name, use_all_sexes = TRUE)
  }
  
  # Observed data
  pts   <- observed_points(test_data, sex, x_var, fixed_var, fixed_val, window)
  n_tot <- attr(pts, "n_total")
  int_scores <- nrow(pts) > 0 && all(abs(pts$y - round(pts$y)) < 1e-8)
  jitter_h   <- if (int_scores) 0.2 else 0
  pts_note <- if (n_tot > nrow(pts)) {
    sprintf("random %s of %s shown", fmt_n(nrow(pts)), fmt_n(n_tot))
  } else {
    sprintf("n = %s", fmt_n(n_tot))
  }
  caption <- sprintf(
    "Band: z = \u00B1%s (centiles %.1f\u2013%.1f) of the fitted %s distribution.\nPoints: %s (%s).",
    k, 100 * pnorm(-k), 100 * pnorm(k), fam, point_desc, pts_note
  )
  
  lbl <- band_label_positions(df, k, y_range)
  
  ggplot(df, aes(x = x)) +
    points_layer(pts, jitter_w, jitter_h, col) +
    geom_ribbon(aes(ymin = pmax(lo, y_range[1]), ymax = pmin(hi, y_range[2])),
                fill = col, alpha = 0.2) +
    geom_line(aes(y = lo),  colour = col, linewidth = 0.5, linetype = "dashed") +
    geom_line(aes(y = hi),  colour = col, linewidth = 0.5, linetype = "dashed") +
    geom_line(aes(y = med), colour = col, linewidth = 1.1) +
    scale_y_continuous(
      sec.axis = dup_axis(breaks = lbl$y, labels = lbl$lab, name = NULL)
    ) +
    coord_cartesian(xlim = xlim, ylim = y_range) +
    labs(
      title    = paste0(test_name, " \u2014 ", sex),
      subtitle = subtitle,
      x        = x_lab,
      y        = score_label(test_name),
      caption  = NULL
    ) +
    theme_pub() +
    theme(
      axis.text.y.right  = element_text(colour = col, face = "bold", size = rel(0.85)),
      axis.ticks.y.right = element_blank()
    )
}

slice_by_age <- function(test_name, long_data, models,
                         sex       = "Male",
                         fixed_age = CROSS_SECTION_AGE,
                         sd_levels = SD_LEVELS,
                         y_range   = NULL) {
  plot_slice(
    test_name, long_data, models, sex,
    x_var      = "Edu",
    x_seq      = seq(min(EDU_RANGE), max(EDU_RANGE), length.out = 200),
    fixed_var  = "age_at_test", fixed_val = fixed_age, window = AGE_WINDOW,
    col        = "#2E6F8E",
    x_lab      = "Education Rank",
    xlim       = range(EDU_RANGE),
    subtitle   = paste0("Score ~ Education (Age = ", fixed_age, ")"),
    point_desc = sprintf("observed baseline, age %s\u2013%s y",
                         fixed_age - AGE_WINDOW, fixed_age + AGE_WINDOW),
    sd_levels  = sd_levels, y_range = y_range,
    jitter_w   = 0.12   # education is discrete
  )
}

slice_by_edu <- function(test_name, long_data, models,
                         sex       = "Male",
                         fixed_edu = CROSS_SECTION_EDU,
                         sd_levels = SD_LEVELS,
                         y_range   = NULL) {
  plot_slice(
    test_name, long_data, models, sex,
    x_var      = "age_at_test",
    x_seq      = seq(AGE_RANGE[1], AGE_RANGE[2], length.out = 200),
    fixed_var  = "Edu", fixed_val = fixed_edu, window = 0,
    col        = "#29AF7F",
    x_lab      = "Age (years)",
    xlim       = AGE_RANGE,
    subtitle   = paste0("Score ~ Age (Edu = ", fixed_edu, ")"),
    point_desc = sprintf("observed baseline, Edu = %s", fixed_edu),
    sd_levels  = sd_levels, y_range = y_range,
    jitter_w   = 0
  )
}

# SAVE OUTPUTS

save_html <- function(p, filename, subdir = NULL) {
  outdir <- if (!is.null(subdir)) file.path(OUTDIR, subdir) else OUTDIR
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  fpath <- file.path(outdir, filename)
  htmlwidgets::saveWidget(as_widget(p), fpath, selfcontained = TRUE)
  sc <- sub("\\.html$", "_files", fpath)
  if (dir.exists(sc)) unlink(sc, recursive = TRUE, force = TRUE)
  message("Saved: ", fpath)
  fpath
}

save_svg <- function(p, filename, subdir = NULL, width = 8, height = 5) {
  outdir <- if (!is.null(subdir)) file.path(OUTDIR, subdir) else OUTDIR
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  fpath <- file.path(outdir, filename)
  ggsave(fpath, plot = p, device = svglite::svglite,
         width = width, height = height, dpi = 300)
  message("Saved: ", fpath)
  fpath
}

save_png <- function(html_file) {
  f_png <- sub("\\.html$", ".png", html_file)
  tryCatch({
    webshot2::webshot(html_file, f_png, vwidth = 800, vheight = 700, zoom = 2)
    message("Saved: ", f_png)
  }, error = function(e) {
    message("Warning: Could not save PNG - ", e$message)
  })
}

# HIGH-RESOLUTION STATIC SURFACE EXPORT
# Plotly 3D surfaces rasterise on export, so use a high-scale PNG (webshot2).
# Render at the SAME canvas size as the layout (SURF_WIDTH x SURF_HEIGHT)
# and raise `scale`:
# fonts keep their intended proportion and resolution goes up. (Rendering on
# a larger canvas instead makes the text proportionally smaller.)
# Defaults give 3200 x 2680 px.
save_surface_image <- function(p, filename, subdir = NULL,
                               width = SURF_WIDTH, height = SURF_HEIGHT,
                               scale = 4) {
  outdir <- if (!is.null(subdir)) file.path(OUTDIR, subdir) else OUTDIR
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  fpath <- file.path(outdir, filename)
  tmp_html <- tempfile(fileext = ".html")
  htmlwidgets::saveWidget(plotly::as_widget(p), tmp_html, selfcontained = TRUE)
  webshot2::webshot(tmp_html, fpath, vwidth = width, vheight = height, zoom = scale)
  unlink(tmp_html)
  unlink(sub("\\.html$", "_files", tmp_html), recursive = TRUE)
  message("Saved: ", fpath)
  fpath
}

# TEST RUNNER
# For one test: load models, then render the 3D surface and both slices for
# each sex. Outputs go to OUTDIR/<test>/.

run_for_test <- function(test) {
  models <- tryCatch(
    load_models_for_test(test),
    error = function(e) { message("[SKIP] ", test, ": ", e$message); NULL }
  )
  if (is.null(models)) return(invisible(FALSE))
  
  z_rng <- SCORE_RANGES[[test]]
  
  for (sx in SEXES) {
    p_surf <- visualise_surface(
      test_name = test, long_data = long_data, models = models,
      age_range = AGE_RANGE, edu_range = EDU_RANGE,
      z_range = z_rng, n_points = N_POINTS,
      sex = sx, sd_levels = SD_LEVELS
    )
    save_html(p_surf, sprintf("%s_%s_surface.html", test, sx), subdir = test)
    save_surface_image(p_surf, sprintf("%s_%s_surface.png", test, sx),
                       subdir = test)
    
    p_age <- slice_by_age(test, long_data, models,
                          sex = sx, fixed_age = CROSS_SECTION_AGE, y_range = z_rng)
    save_svg(p_age, sprintf("%s_%s_slice_age%s.svg", test, sx, CROSS_SECTION_AGE),
             subdir = test)
    
    p_edu <- slice_by_edu(test, long_data, models,
                          sex = sx, fixed_edu = CROSS_SECTION_EDU, y_range = z_rng)
    save_svg(p_edu, sprintf("%s_%s_slice_edu%s.svg", test, sx, CROSS_SECTION_EDU),
             subdir = test)
  }
  
  invisible(TRUE)
}

# CLI ENTRY POINT
#   Rscript vis_batch.R              # all tests
#   Rscript vis_batch.R --test BL    # single test

args      <- commandArgs(trailingOnly = TRUE)
tests_all <- names(OPTIMAL_DISTRIBUTIONS)

if (length(args) >= 2 && args[1] == "--test") {
  test <- args[2]
  if (!test %in% tests_all) stop("Unknown test: ", test)
  ok <- run_for_test(test)
  message("Done: ", test, " (ok=", ok, ")")
} else {
  n_ok <- 0L; n_skip <- 0L
  for (t in tests_all) {
    ok <- run_for_test(t)
    if (isTRUE(ok)) n_ok <- n_ok + 1L else n_skip <- n_skip + 1L
  }
  message("----------- Summary -----------")
  message("Rendered: ", n_ok, " test(s)")
  message("Skipped : ", n_skip, " test(s)")
  message("Output  : ", normalizePath(OUTDIR, mustWork = FALSE))
}