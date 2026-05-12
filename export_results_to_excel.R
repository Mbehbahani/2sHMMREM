# =============================================================================
# Export HMMREM Results to Excel
# This script runs the analysis and exports all results to an Excel file
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Setup
# -----------------------------------------------------------------------------

# Skip package loading if already loaded by runner script
if (!exists("PACKAGES_LOADED_BY_RUNNER") || !PACKAGES_LOADED_BY_RUNNER) {
  # Check and install required packages
  required_packages <- c("dplyr", "ggplot2", 
                         "plotly", "momentuHMM", "remulate", "tidyr", 
                         "knitr", "openxlsx")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg == "openxlsx") {
        message("Installing openxlsx package for Excel export...")
        install.packages("openxlsx")
      } else {
        stop(paste("Package", pkg, "is required but not installed."))
      }
    }
    library(pkg, character.only = TRUE)
  }
}

load_analysis_bundle <- function(bundle_path) {
  bundle <- readRDS(bundle_path)

  if (!is.list(bundle)) {
    stop("Analysis bundle must be a named list.")
  }

  required_names <- c(
    "all_results", "all_hmm_agg", "all_diagnostics", "all_rem_agg",
    "scenarios", "run_configs", "R", "n_actors", "states_to_fit",
    "tau1", "tau2", "RUN_HMM_REM"
  )
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

# -----------------------------------------------------------------------------
# 2. Run the main analysis (source the scenarios script) if not already done
# -----------------------------------------------------------------------------

bundle_path <- "hmmrem_analysis_bundle.rds"

if (!file.exists(bundle_path)) {
  cat("Analysis bundle not found. Running HMMREM Scenario Analysis...\n")
  cat("This may take several minutes.\n\n")

  source("HMMREM_scenarios.R")

  if (!file.exists(bundle_path)) {
    stop(
      "Analysis completed without creating ",
      bundle_path,
      ". Please check ",
      "HMMREM_scenarios.R",
      " for bundle-saving errors."
    )
  }
}

cat("Loading analysis bundle from:", normalizePath(bundle_path, winslash = "/"), "\n\n")
load_analysis_bundle(bundle_path)

# Check what results are available
# HMM 2-state is ALWAYS fitted (even in quick mode), but require at least one non-empty aggregate
hmm_available <- FALSE
if (exists("all_hmm_agg") && length(all_hmm_agg) > 0) {
  valid_hmm_agg <- Filter(function(x) !is.null(x) && nrow(x) > 0, all_hmm_agg)
  hmm_available <- length(valid_hmm_agg) > 0
}
if (!hmm_available) {
  cat("  DEBUG: exists('all_hmm_agg') =", exists("all_hmm_agg"), "\n")
  if (exists("all_hmm_agg")) {
    cat("  DEBUG: length(all_hmm_agg)  =", length(all_hmm_agg), "\n")
    cat("  DEBUG: names               =", paste(head(names(all_hmm_agg), 5), collapse = ", "), "\n")
    if (length(all_hmm_agg) > 0) {
      cat("  DEBUG: class of first entry =", class(all_hmm_agg[[1]]), "\n")
      if (is.data.frame(all_hmm_agg[[1]])) cat("  DEBUG: nrow of first entry  =", nrow(all_hmm_agg[[1]]), "\n")
    }
  }
}

rem_available <- FALSE
if (exists("all_rem_agg") && length(all_rem_agg) > 0) {
  valid_rem_agg <- Filter(function(x) !is.null(x) && length(x) > 0, all_rem_agg)
  rem_available <- length(valid_rem_agg) > 0
}

extract_rem_uncertainty <- function(all_rem_agg, run_configs) {
  if (!exists("all_rem_agg") || length(all_rem_agg) == 0) return(data.frame())

  rows <- lapply(names(all_rem_agg), function(config_name) {
    config_res <- all_rem_agg[[config_name]]
    if (is.null(config_res) || length(config_res) == 0) return(NULL)

    model_rows <- lapply(names(config_res), function(model_name) {
      model_res <- config_res[[model_name]]
      coef_df <- model_res$coefficients
      if (is.null(coef_df) || nrow(coef_df) == 0) return(NULL)

      data.frame(
        Config = config_name,
        Scenario = run_configs[[config_name]]$scenario_name,
        T = run_configs[[config_name]]$n_events,
        N = run_configs[[config_name]]$n_actors,
        Model = model_name,
        Variable = coef_df$Variable,
        Estimate_mean = coef_df$Estimate_mean,
        Empirical_SD = coef_df$Estimate_sd,
        Reported_SE = coef_df$SE_mean,
        Coverage_95 = coef_df$Coverage_95,
        SE_to_SD_ratio = coef_df$SE_to_SD_ratio,
        SD_to_SE_ratio = coef_df$SD_to_SE_ratio,
        sig_rate = coef_df$sig_rate,
        BIC_mean = model_res$bic_mean,
        BIC_sd = model_res$bic_sd
      )
    })

    model_rows <- Filter(Negate(is.null), model_rows)
    if (length(model_rows) == 0) return(NULL)
    do.call(rbind, model_rows)
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}
# -----------------------------------------------------------------------------
# 3. Prepare data for Excel export
# -----------------------------------------------------------------------------

cat("\n\nPreparing Excel export...\n")

# Create workbook
wb <- createWorkbook()

# -----------------------------------------------------------------------------
# Sheet 1: Analysis Summary
# -----------------------------------------------------------------------------

addWorksheet(wb, "Summary")

summary_info <- data.frame(
  Parameter = c("Analysis Date", "Number of Replications", 
                "Actor Counts", "States Fitted", "Confidence Threshold 1 (tau1)",
                "Confidence Threshold 2 (tau2)",
                "Scenarios", "Configurations", "REM Models"),
  Value = c(as.character(Sys.time()), R, paste(n_actors, collapse = ", "), 
            if (exists("RUN_HMM_REM") && RUN_HMM_REM) paste(states_to_fit, collapse = ", ") else "2", tau1, tau2,
            paste(names(scenarios), collapse = ", "),
            paste(length(run_configs), "total"),
            if (rem_available) paste(unique(unlist(lapply(Filter(function(x) !is.null(x) && length(x) > 0, all_rem_agg), names))), collapse = ", ") else "None")
)

writeData(wb, "Summary", summary_info, startRow = 1, startCol = 1)

# Add scenario parameters
addWorksheet(wb, "Scenario_Parameters")

# Shared configuration lookup (used across multiple tabs)
config_lookup <- data.frame(
  Config = names(run_configs),
  Scenario = sapply(run_configs, function(x) x$scenario_name),
  T = sapply(run_configs, function(x) x$n_events),
  N = sapply(run_configs, function(x) x$n_actors),
  stringsAsFactors = FALSE
)

# Use run_configs for scenario parameters (scenario x T x N combinations)
scenario_params <- config_lookup %>%
  mutate(
  transition_stay_state1 = sapply(run_configs, function(x) x$scenario$transition_probs[1,1]),
  transition_stay_state2 = sapply(run_configs, function(x) x$scenario$transition_probs[2,2]),
  sep_factor = sapply(run_configs, function(x) x$scenario$sep_factor),
  dwell_state1 = sapply(run_configs, function(x) round(x$scenario$dwell_time[1], 2)),
  dwell_state2 = sapply(run_configs, function(x) round(x$scenario$dwell_time[2], 2)),
  baseline1 = sapply(run_configs, function(x) x$scenario$baseline1),
  baseline2 = sapply(run_configs, function(x) x$scenario$baseline2),
  outdegree1 = sapply(run_configs, function(x) x$scenario$outdegree1),
  outdegree2 = sapply(run_configs, function(x) x$scenario$outdegree2),
  inertia1 = sapply(run_configs, function(x) x$scenario$inertia1),
  inertia2 = sapply(run_configs, function(x) x$scenario$inertia2)
)

writeData(wb, "Scenario_Parameters", scenario_params, startRow = 1, startCol = 1)

# -----------------------------------------------------------------------------
# Sheet 3: Diagnostics (rho_hat, dwell_time, sep_factor)
# -----------------------------------------------------------------------------

addWorksheet(wb, "Diagnostics")

if (exists("all_diagnostics") && length(all_diagnostics) > 0) {
  diag_df <- do.call(rbind, lapply(all_diagnostics, function(d) {
    data.frame(
      Config = d$config,
      Scenario = d$scenario,
      T = d$T,
      N = d$n_actors,
      sep_factor = d$sep_factor,
      Dwell_State1 = round(d$dwell_time[1], 2),
      Dwell_State2 = round(d$dwell_time[2], 2),
      rho_hat_mean = round(d$rho_hat_mean, 4),
      rho_hat_sd = round(d$rho_hat_sd, 4)
    )
  }))
  writeData(wb, "Diagnostics", diag_df, startRow = 1, startCol = 1)
}

# -----------------------------------------------------------------------------
# Sheet 4: HMM Results (always available - at least 2-state is fitted)
# -----------------------------------------------------------------------------

if (hmm_available) {

addWorksheet(wb, "HMM_Results")

# Prepare hmm_combined data (needed for multiple sheets)
valid_hmm_agg <- Filter(function(x) !is.null(x) && nrow(x) > 0, all_hmm_agg)
all_hmm_results <- if (length(valid_hmm_agg) > 0) do.call(rbind, valid_hmm_agg) else data.frame()
hmm_combined <- all_hmm_results %>%
  mutate(Config = as.character(Scenario)) %>%
  select(-Scenario) %>%
  left_join(config_lookup, by = "Config") %>%
  select(Config, Scenario, T, N, everything())

writeData(wb, "HMM_Results", hmm_combined, startRow = 1, startCol = 1) 

# Add conditional formatting for best BIC
if (nrow(hmm_combined) > 0 && "BIC_mean" %in% colnames(hmm_combined)) {
  conditionalFormatting(wb, "HMM_Results", 
                        cols = which(colnames(hmm_combined) == "BIC_mean"),
                        rows = 2:(nrow(hmm_combined) + 1),
                        type = "colourScale",
                        style = c("#63BE7B", "#FFEB84", "#F8696B"))
}

# -----------------------------------------------------------------------------
# Sheet 5: Confidence Metrics (always available - 2-state HMM is always fitted)
# -----------------------------------------------------------------------------

addWorksheet(wb, "Confidence_Metrics")

confidence_df <- hmm_combined %>%
  filter(States == 2) %>%
  select(Config, Scenario, T, N, Confidence_mean, Confidence_sd, 
         Uncertainty_tau1_mean, Uncertainty_tau1_sd,
         Uncertainty_tau2_mean, Uncertainty_tau2_sd,
         Accuracy_mean, Accuracy_sd)

writeData(wb, "Confidence_Metrics", confidence_df, startRow = 1, startCol = 1)

# -----------------------------------------------------------------------------
# Sheet 6: State Overlap Diagnostics (2-state only)
# -----------------------------------------------------------------------------

addWorksheet(wb, "State_Overlap")

overlap_cols <- c("OVL_hist_mean", "OVL_hist_sd", "OVL_kde_mean", "OVL_kde_sd",
                  "W1_mean", "W1_sd", "symKL_param_mean", "symKL_param_sd",
                  "OVL_param_mean", "OVL_param_sd", "JS_param_mean", "JS_param_sd",
                  "n_valid_OVL_hist", "n_valid_OVL_kde", "n_valid_W1",
                  "n_valid_symKL_param", "n_valid_OVL_param", "n_valid_JS_param")
# Only include overlap columns that exist in the data
existing_overlap_cols <- intersect(overlap_cols, colnames(hmm_combined))

if (length(existing_overlap_cols) > 0) {
  overlap_df <- hmm_combined %>%
    filter(States == 2) %>%
    select(Config, Scenario, T, N, all_of(existing_overlap_cols))
  writeData(wb, "State_Overlap", overlap_df, startRow = 1, startCol = 1)
  
  # Add conditional formatting for OVL columns (green=low overlap, red=high)
  for (col_name in c("OVL_hist_mean", "OVL_kde_mean")) {
    if (col_name %in% colnames(overlap_df)) {
      col_idx <- which(colnames(overlap_df) == col_name)
      if (nrow(overlap_df) > 0) {
        conditionalFormatting(wb, "State_Overlap",
                              cols = col_idx,
                              rows = 2:(nrow(overlap_df) + 1),
                              type = "colourScale",
                              style = c("#63BE7B", "#FFEB84", "#F8696B"))
      }
    }
  }
} else {
  writeData(wb, "State_Overlap",
            data.frame(Note = "No overlap data available. Re-run analysis with updated HMMREM_scenarios.R."),
            startRow = 1, startCol = 1)
}

} else {
  cat("  No HMM results available. Sheets 4-6 (HMM_Results, Confidence_Metrics, State_Overlap) skipped.\n")
}

if (rem_available) {
  addWorksheet(wb, "REM_Uncertainty")
  rem_uncertainty_df <- extract_rem_uncertainty(all_rem_agg, run_configs)

  if (nrow(rem_uncertainty_df) > 0) {
    writeData(wb, "REM_Uncertainty", rem_uncertainty_df, startRow = 1, startCol = 1)

    if ("SE_to_SD_ratio" %in% colnames(rem_uncertainty_df)) {
      ratio_col <- which(colnames(rem_uncertainty_df) == "SE_to_SD_ratio")
      conditionalFormatting(
        wb, "REM_Uncertainty",
        cols = ratio_col,
        rows = 2:(nrow(rem_uncertainty_df) + 1),
        type = "colourScale",
        style = c("#F8696B", "#FFEB84", "#63BE7B")
      )
    }
  } else {
    writeData(
      wb, "REM_Uncertainty",
      data.frame(Note = "REM results object exists, but no uncertainty summary rows were available."),
      startRow = 1, startCol = 1
    )
  }
} else {
  cat("  No REM uncertainty results available. REM_Uncertainty sheet skipped.\n")
}

# -----------------------------------------------------------------------------
# Save Excel file
# -----------------------------------------------------------------------------

excel_filename <- "HMMREM_results.xlsx"
saveWorkbook(wb, excel_filename, overwrite = TRUE)

cat("\n======================================================================\n")
cat("Excel file saved successfully!\n")
cat("Filename:", excel_filename, "\n")
cat("Location:", getwd(), "\n")
cat("======================================================================\n\n")

cat("Sheets in the Excel file:\n")
cat("  1. Summary - Analysis parameters and settings (tau1=0.9, tau2=0.8)\n")
cat("  2. Scenario_Parameters - True parameters for each config\n")
cat("  3. Diagnostics - rho_hat, dwell_time, sep_factor per config\n")

if (hmm_available) {
  cat("  4. HMM_Results - BIC, accuracy, confidence for HMM models\n")
  cat("  5. Confidence_Metrics - State assignment confidence with tau1 & tau2 (2-state model)\n")
  cat("  6. State_Overlap - State overlap diagnostics (2-state model)\n")
}

if (rem_available) {
  cat("  7. REM_Uncertainty - Empirical SD vs reported SE for second-stage REM estimates\n")
}

if (!hmm_available) {
  cat("  (No HMM results: sheets 4-6 skipped)\n")
}
if (!rem_available) {
  cat("  (No REM uncertainty results: REM sheet skipped)\n")
}
