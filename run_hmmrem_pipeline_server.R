# =============================================================================
# Server Runner: HMMREM full pipeline
# - Checks/installs required libraries
# - Runs HMMREM_scenarios.R
# - Runs export_results_to_excel.R
# =============================================================================

# Use Posit Public Package Manager for pre-compiled binaries (Ubuntu 24.04 compatible)
options(repos = c(
  RSPM = "https://packagemanager.posit.co/cran/__linux__/noble/latest",
  CRAN = "https://cloud.r-project.org"
))

cat("\n============================================================\n")
cat("HMMREM SERVER PIPELINE START\n")
cat("Working directory:", normalizePath(getwd(), winslash = "/"), "\n")
cat("R version:", R.version.string, "\n")
cat("Start time:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

required_packages <- c(
  "remstats", "remstimate", "dplyr", "ggplot2", "plotly",
  "momentuHMM", "remulate", "tidyr", "knitr", "openxlsx", "remify"
)

install_if_missing <- function(pkgs) {
  cat("Checking package installation status...\n")
  installed <- rownames(installed.packages())
  missing <- setdiff(pkgs, installed)
  
  if (length(missing) > 0) {
    cat("Missing packages:", paste(missing, collapse = ", "), "\n\n")
    
    for (pkg in missing) {
      cat("Installing", pkg, "... ")
      result <- tryCatch({
        # Force binary installation first (much faster, no compilation)
        install.packages(pkg, dependencies = TRUE, quiet = TRUE, 
                        INSTALL_opts = "--no-test-load")
        if (requireNamespace(pkg, quietly = TRUE)) {
          cat("OK\n")
          TRUE
        } else {
          cat("FAILED\n")
          FALSE
        }
      }, error = function(e) {
        cat("ERROR:", conditionMessage(e), "\n")
        FALSE
      })
      
      if (!result) {
        stop(paste("Failed to install required package:", pkg, 
                  "\nTry running manually: install.packages('", pkg, "')"))
      }
    }
    cat("\nAll missing packages installed successfully.\n")
  } else {
    cat("All required packages are already installed.\n")
  }
  
  cat("\nLoading packages...\n")
  for (pkg in pkgs) {
    cat("  Loading", pkg, "... ")
    # Only load if not already loaded
    if (!isNamespaceLoaded(pkg)) {
      if (!suppressPackageStartupMessages(library(pkg, character.only = TRUE, 
                                                  logical.return = TRUE))) {
        cat("FAILED\n")
        stop(paste("Failed to load package:", pkg))
      }
      cat("OK\n")
    } else {
      cat("ALREADY LOADED\n")
    }
  }
  cat("All packages loaded successfully.\n")
  
  # Set a flag to tell sourced scripts that packages are already loaded
  assign("PACKAGES_LOADED_BY_RUNNER", TRUE, envir = .GlobalEnv)
}

run_script <- function(script_path, label) {
  if (!file.exists(script_path)) {
    stop(paste0("Required script not found: ", script_path))
  }
  
  cat("\n------------------------------------------------------------\n")
  cat("Running:", label, "\n")
  cat("Script:", script_path, "\n")
  cat("------------------------------------------------------------\n")
  
  t0 <- Sys.time()
  source(script_path, echo = FALSE)
  t1 <- Sys.time()
  
  cat("Completed:", label, "in", round(as.numeric(difftime(t1, t0, units = "secs")), 2), "seconds\n")
}

tryCatch({
  install_if_missing(required_packages)
  
  run_script("HMMREM_scenarios.R", "Scenario simulation + model fitting")
  run_script("export_results_to_excel.R", "Excel export")
  
  cat("\n============================================================\n")
  cat("HMMREM SERVER PIPELINE FINISHED SUCCESSFULLY\n")
  cat("Output file expected: HMMREM_results.xlsx\n")
  cat("End time:", as.character(Sys.time()), "\n")
  cat("============================================================\n\n")
  
}, error = function(e) {
  cat("\n============================================================\n")
  cat("PIPELINE FAILED\n")
  cat("Error:", conditionMessage(e), "\n")
  cat("End time:", as.character(Sys.time()), "\n")
  cat("============================================================\n\n")
  quit(save = "no", status = 1)
})
