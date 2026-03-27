# =============================================================================
# Rebuild plots and Excel export from saved HMMREM outputs
# - Recreates per-configuration density plots from density_data/*.rds
# - Recreates comparison plots when HMM summaries are available
# - Reuses export_results_to_excel.R without re-running HMMREM_scenarios.R
#
# Usage:
#   Rscript create_plots_and_excel_from_saved_results.R
#   Rscript create_plots_and_excel_from_saved_results.R path/to/hmmrem_analysis_bundle.rds
#
# Notes:
# - density_data/*.rds is enough for the density plots only.
# - The Excel workbook and HMM comparison plots also need all_hmm_agg and
#   all_diagnostics, either already in the workspace or stored in a bundle .rds.
# =============================================================================

required_packages <- c("dplyr", "ggplot2", "tidyr", "openxlsx")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed."))
  }
  library(pkg, character.only = TRUE)
}

build_hmmrem_context <- function() {
  base_emissions <- list(
    baseline1 = -8, baseline2 = -6.5,
    outdegree1 = 0.2, outdegree2 = 0.5,
    inertia1 = 0.2, inertia2 = 0.4
  )

  compute_scaled_emissions <- function(base, sep_factor) {
    params <- c("baseline", "outdegree", "inertia")
    result <- list()

    for (p in params) {
      p1 <- base[[paste0(p, "1")]]
      p2 <- base[[paste0(p, "2")]]
      mid <- (p1 + p2) / 2
      half_diff <- (p2 - p1) / 2
      result[[paste0(p, "1")]] <- mid - sep_factor * half_diff
      result[[paste0(p, "2")]] <- mid + sep_factor * half_diff
    }

    result
  }

  sep_factor_hard <- 0.75
  hard_emissions <- compute_scaled_emissions(base_emissions, sep_factor_hard)

  scenarios <- list(
    Easy = list(
      name = "Easy",
      transition_probs = matrix(c(0.99, 0.01, 0.02, 0.98), nrow = 2, byrow = TRUE),
      baseline1 = base_emissions$baseline1,
      baseline2 = base_emissions$baseline2,
      outdegree1 = base_emissions$outdegree1,
      outdegree2 = base_emissions$outdegree2,
      inertia1 = base_emissions$inertia1,
      inertia2 = base_emissions$inertia2,
      sep_factor = 1.0,
      T_values = c(3000, 6000)
    ),
    Medium = list(
      name = "Medium",
      transition_probs = matrix(c(0.9, 0.1, 0.15, 0.85), nrow = 2, byrow = TRUE),
      baseline1 = base_emissions$baseline1,
      baseline2 = base_emissions$baseline2,
      outdegree1 = base_emissions$outdegree1,
      outdegree2 = base_emissions$outdegree2,
      inertia1 = base_emissions$inertia1,
      inertia2 = base_emissions$inertia2,
      sep_factor = 1.0,
      T_values = c(3000, 6000)
    ),
    Hard = list(
      name = "Hard",
      transition_probs = matrix(c(0.9, 0.1, 0.15, 0.85), nrow = 2, byrow = TRUE),
      baseline1 = hard_emissions$baseline1,
      baseline2 = hard_emissions$baseline2,
      outdegree1 = hard_emissions$outdegree1,
      outdegree2 = hard_emissions$outdegree2,
      inertia1 = hard_emissions$inertia1,
      inertia2 = hard_emissions$inertia2,
      sep_factor = sep_factor_hard,
      T_values = c(3000, 6000)
    ),
    ExtremeHard = list(
      name = "ExtremeHard",
      transition_probs = matrix(c(0.75, 0.25, 0.30, 0.70), nrow = 2, byrow = TRUE),
      baseline1 = hard_emissions$baseline1,
      baseline2 = hard_emissions$baseline2,
      outdegree1 = hard_emissions$outdegree1,
      outdegree2 = hard_emissions$outdegree2,
      inertia1 = hard_emissions$inertia1,
      inertia2 = hard_emissions$inertia2,
      sep_factor = sep_factor_hard,
      T_values = c(3000, 6000)
    )
  )

  for (sc_name in names(scenarios)) {
    persistence <- diag(scenarios[[sc_name]]$transition_probs)
    scenarios[[sc_name]]$dwell_time <- 1 / (1 - persistence)
  }

  R <- 100
  n_actors <- c(10, 15, 20)
  states_to_fit <- 1:4
  tau1 <- 0.9
  tau2 <- 0.8
  RUN_HMM_REM <- TRUE

  run_configs <- list()
  for (sc_name in names(scenarios)) {
    sc <- scenarios[[sc_name]]
    for (T_val in sc$T_values) {
      for (n_val in n_actors) {
        config_name <- paste0(sc_name, "_T", T_val, "_N", n_val)
        run_configs[[config_name]] <- list(
          scenario_name = sc_name,
          config_name = config_name,
          n_events = T_val,
          n_actors = n_val,
          scenario = sc
        )
      }
    }
  }

  list(
    scenarios = scenarios,
    run_configs = run_configs,
    R = R,
    n_actors = n_actors,
    states_to_fit = states_to_fit,
    tau1 = tau1,
    tau2 = tau2,
    RUN_HMM_REM = RUN_HMM_REM
  )
}

create_density_plot <- function(density_data) {
  stopifnot("Data frame must contain 'x' column" = "x" %in% colnames(density_data))
  stopifnot("Data frame must contain 'y' column" = "y" %in% colnames(density_data))
  stopifnot("Data frame must contain 'State' column" = "State" %in% colnames(density_data))
  stopifnot("Data frame cannot be empty" = nrow(density_data) > 0)
  stopifnot("'x' must be numeric" = is.numeric(density_data$x))
  stopifnot("'y' must be numeric" = is.numeric(density_data$y))

  if (!"replication" %in% colnames(density_data)) {
    warning("Missing 'replication' column; falling back to raw density lines.")
    return(
      ggplot(density_data, aes(x = x, y = y, color = State)) +
        geom_line(linewidth = 1) +
        labs(x = "Inter-event time", y = "Density") +
        theme_minimal() +
        theme(legend.position = "bottom")
    )
  }

  common_x <- seq(0, max(density_data$x, na.rm = TRUE), length.out = 512)
  density_groups <- split(
    density_data,
    list(density_data$State, density_data$replication),
    drop = TRUE
  )

  interpolated_density <- lapply(density_groups, function(df) {
    df <- df[order(df$x), c("x", "y", "State", "replication")]
    df <- df[!duplicated(df$x), ]
    interpolated_y <- stats::approx(
      x = df$x,
      y = df$y,
      xout = common_x,
      yleft = 0,
      yright = 0,
      rule = 1,
      ties = mean
    )$y

    data.frame(
      x = common_x,
      y = interpolated_y,
      State = df$State[1],
      replication = df$replication[1]
    )
  })

  summary_density <- dplyr::bind_rows(interpolated_density) |>
    dplyr::group_by(State, x) |>
    dplyr::summarise(
      density_mean = mean(y, na.rm = TRUE),
      density_min = min(y, na.rm = TRUE),
      density_max = max(y, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot(summary_density, aes(x = x, color = State, fill = State)) +
    geom_ribbon(
      aes(ymin = density_min, ymax = density_max),
      alpha = 0.35,
      color = NA,
      show.legend = FALSE
    ) +
    geom_line(aes(y = density_mean), linewidth = 0.7) +
    labs(x = "Inter-event time", y = "Density") +
    scale_x_continuous(limits = c(0, 60), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 60)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.border = element_rect(colour = "#00000080", fill = NA, linewidth = 1.5),
      text = element_text(size = 16),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_blank()
    )
}

create_bic_comparison_plot <- function(all_hmm_results, R_val, config_levels) {
  combined <- do.call(rbind, all_hmm_results)
  combined$Scenario <- factor(combined$Scenario, levels = config_levels)
  combined$States <- factor(combined$States)

  ggplot(combined, aes(x = States, y = BIC_mean, fill = Scenario)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(
      aes(ymin = BIC_mean - BIC_sd, ymax = BIC_mean + BIC_sd),
      position = position_dodge(width = 0.8),
      width = 0.2
    ) +
    labs(
      title = "HMM Model Comparison: BIC by Number of States",
      subtitle = paste0("Based on ", R_val, " replications per scenario"),
      x = "Number of Hidden States",
      y = "BIC (mean +/- SD)",
      fill = "Scenario Difficulty"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.position = "bottom"
    )
}

create_confidence_plot <- function(all_hmm_results, tau1_val, tau2_val, config_levels) {
  combined <- do.call(rbind, all_hmm_results)
  combined$Scenario <- factor(combined$Scenario, levels = config_levels)
  combined_2state <- combined %>% filter(States == 2)

  metric_tau1 <- paste0("% Uncertain (tau=", tau1_val, ")")
  metric_tau2 <- paste0("% Uncertain (tau=", tau2_val, ")")

  plot_data <- combined_2state %>%
    select(Scenario, Confidence_mean, Uncertainty_tau1_mean, Uncertainty_tau2_mean) %>%
    mutate(Confidence_pct = Confidence_mean * 100) %>%
    pivot_longer(
      cols = c(Confidence_pct, Uncertainty_tau1_mean, Uncertainty_tau2_mean),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(Metric = case_when(
      Metric == "Confidence_pct" ~ "Avg Confidence",
      Metric == "Uncertainty_tau1_mean" ~ metric_tau1,
      Metric == "Uncertainty_tau2_mean" ~ metric_tau2
    ))

  color_map <- c(
    "Avg Confidence" = "#3498DB",
    metric_tau1 = "#E74C3C",
    metric_tau2 = "#F39C12"
  )

  ggplot(plot_data, aes(x = Scenario, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = color_map) +
    labs(
      title = "State Classification Confidence (2-State HMM)",
      subtitle = paste0("Uncertainty thresholds: tau1=", tau1_val, ", tau2=", tau2_val),
      x = "Scenario Difficulty",
      y = "Percentage (%)",
      fill = "Metric"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "bottom"
    )
}

create_overlap_plot <- function(all_hmm_results) {
  combined <- do.call(rbind, all_hmm_results)
  combined_2state <- combined %>% filter(States == 2)

  if (nrow(combined_2state) == 0 ||
      (all(is.na(combined_2state$OVL_hist_mean)) && all(is.na(combined_2state$OVL_kde_mean)))) {
    return(NULL)
  }

  plot_data <- combined_2state %>%
    select(Scenario, OVL_hist_mean, OVL_kde_mean) %>%
    pivot_longer(
      cols = c(OVL_hist_mean, OVL_kde_mean),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(Metric = case_when(
      Metric == "OVL_hist_mean" ~ "Histogram OVL",
      Metric == "OVL_kde_mean" ~ "KDE OVL"
    ))

  ggplot(plot_data, aes(x = Scenario, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "State Overlap Coefficients (2-State HMM)",
      subtitle = "Higher overlap = harder to distinguish states",
      x = "Configuration",
      y = "Overlap Coefficient",
      fill = "Method"
    ) +
    scale_fill_manual(values = c("Histogram OVL" = "#E74C3C", "KDE OVL" = "#3498DB")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

load_analysis_bundle <- function(bundle_path) {
  bundle <- readRDS(bundle_path)

  if (!is.list(bundle)) {
    stop("Analysis bundle must be a named list.")
  }

  required_names <- c("all_hmm_agg", "all_diagnostics")
  missing_names <- setdiff(required_names, names(bundle))
  if (length(missing_names) > 0) {
    stop(
      "Analysis bundle is missing required objects: ",
      paste(missing_names, collapse = ", ")
    )
  }

  for (obj_name in names(bundle)) {
    assign(obj_name, bundle[[obj_name]], envir = .GlobalEnv)
  }
}

has_valid_hmm_agg <- function() {
  exists("all_hmm_agg", envir = .GlobalEnv) &&
    length(all_hmm_agg) > 0 &&
    length(Filter(function(x) !is.null(x) && nrow(x) > 0, all_hmm_agg)) > 0
}

has_valid_diagnostics <- function() {
  exists("all_diagnostics", envir = .GlobalEnv) &&
    length(all_diagnostics) > 0
}

write_density_plots <- function(density_dir = "density_data", plot_dir = file.path("plots", "density")) {
  density_files <- list.files(
    density_dir,
    pattern = "^density_data_.*\\.rds$",
    full.names = TRUE
  )

  if (length(density_files) == 0) {
    stop("No density_data_*.rds files were found in ", density_dir)
  }

  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  saved_paths <- character(length(density_files))

  for (i in seq_along(density_files)) {
    density_path <- density_files[[i]]
    density_df <- readRDS(density_path)
    config_name <- unique(density_df$config_name)

    if (length(config_name) != 1 || is.na(config_name)) {
      stop("Could not infer a single config_name from ", density_path)
    }

    p_density <- create_density_plot(density_df)
    output_path <- file.path(plot_dir, paste0("density_", config_name, ".svg"))
    ggsave(output_path, p_density, device = "svg", width = 8, height = 5, bg = "white")
    saved_paths[[i]] <- output_path
  }

  saved_paths
}

write_comparison_plots <- function() {
  if (!has_valid_hmm_agg()) {
    return(invisible(NULL))
  }

  full_config_names <- names(run_configs)[vapply(
    run_configs,
    function(cfg) identical(cfg$n_actors, 15) && identical(cfg$n_events, 3000),
    logical(1)
  )]
  full_hmm_agg <- all_hmm_agg[full_config_names[full_config_names %in% names(all_hmm_agg)]]

  if (length(full_hmm_agg) == 0) {
    message("Skipping comparison plots: no full-analysis HMM summaries were found.")
    return(invisible(NULL))
  }

  p_bic <- create_bic_comparison_plot(full_hmm_agg, R, full_config_names)
  ggsave("HMMREM_scenario_comparison.png", p_bic, width = 10, height = 6, dpi = 300, bg = "white")

  p_conf <- create_confidence_plot(full_hmm_agg, tau1, tau2, full_config_names)
  ggsave("HMMREM_confidence_comparison.png", p_conf, width = 8, height = 6, dpi = 300, bg = "white")

  p_ovl <- create_overlap_plot(full_hmm_agg)
  if (!is.null(p_ovl)) {
    ggsave("HMMREM_overlap_comparison.png", p_ovl, width = 10, height = 6, dpi = 300, bg = "white")
  }
}

args <- commandArgs(trailingOnly = TRUE)
bundle_path <- if (length(args) >= 1) args[[1]] else "hmmrem_analysis_bundle.rds"

context <- build_hmmrem_context()
for (obj_name in names(context)) {
  assign(obj_name, context[[obj_name]], envir = .GlobalEnv)
}

if ((!has_valid_hmm_agg() || !has_valid_diagnostics()) && file.exists(bundle_path)) {
  message("Loading analysis bundle from: ", normalizePath(bundle_path, winslash = "/"))
  load_analysis_bundle(bundle_path)
}

message("Rebuilding density plots from density_data/*.rds ...")
density_plot_paths <- write_density_plots()
message("Saved ", length(density_plot_paths), " density plot(s).")

if (has_valid_hmm_agg()) {
  message("Rebuilding HMM comparison plots ...")
  write_comparison_plots()
}

if (!has_valid_hmm_agg() || !has_valid_diagnostics()) {
  stop(
    paste(
      "Density plots were created, but the exact Excel export cannot be rebuilt",
      "from density_data alone. Provide all_hmm_agg and all_diagnostics in the",
      "current R session or pass a bundle .rds as the first argument."
    )
  )
}

if (!exists("all_results", envir = .GlobalEnv)) {
  assign("all_results", list(), envir = .GlobalEnv)
}

message("Rebuilding Excel workbook with export_results_to_excel.R ...")
source("export_results_to_excel.R", echo = FALSE)
