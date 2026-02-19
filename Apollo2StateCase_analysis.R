# =============================================================================
# Apollo 13 HMMREM Analysis with Comprehensive Metrics
# Based on Apollo2StateCase.Rmd with added metrics from HMMREM_scenarios.R
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Setup and Package Loading
# -----------------------------------------------------------------------------

required_packages <- c("remify", "remstats", "remstimate", "dplyr",  
                       "zoo", "ggplot2", "plotly", "fastDummies", 
                       "abind", "momentuHMM", "knitr", "tidyr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# -----------------------------------------------------------------------------
# 2. Define Analysis Parameters
# -----------------------------------------------------------------------------

# Confidence thresholds for uncertainty metrics
tau1 <- 0.9
tau2 <- 0.8

# Number of states to test
states_to_fit <- 1:4

# Train/Test split
N_train <- 4998
N_test <- 500

# -----------------------------------------------------------------------------
# 3. Helper Functions (from HMMREM_scenarios.R)
# -----------------------------------------------------------------------------

# Compute confidence metrics for HMM state assignments
compute_confidence_metrics <- function(hmm_fit, tau1 = 0.9, tau2 = 0.8) {
  tryCatch({
    state_probs <- momentuHMM::stateProbs(hmm_fit)
    max_probs <- apply(state_probs, 1, max)
    avg_confidence <- mean(max_probs, na.rm = TRUE)
    pct_uncertain_tau1 <- mean(max_probs < tau1, na.rm = TRUE) * 100
    pct_uncertain_tau2 <- mean(max_probs < tau2, na.rm = TRUE) * 100
    return(list(
      avg_confidence = avg_confidence, 
      pct_uncertain_tau1 = pct_uncertain_tau1,
      pct_uncertain_tau2 = pct_uncertain_tau2,
      state_probs = state_probs,
      max_probs = max_probs
    ))
  }, error = function(e) {
    return(list(avg_confidence = NA, pct_uncertain_tau1 = NA, pct_uncertain_tau2 = NA,
                state_probs = NULL, max_probs = NULL))
  })
}

# Fit HMM safely with error handling
fit_hmm_safe <- function(hmm_data, n_states) {
  tryCatch({
    dist <- list(step = "exp")
    if (n_states == 1) {
      Par0 <- list(step = 0.3)
    } else {
      Par0 <- list(step = seq(0.1, 0.8, length.out = n_states))
    }
    HMMfit <- fitHMM(data = hmm_data, nbStates = n_states, dist = dist, 
                     Par0 = Par0, formula = ~1)
    return(HMMfit)
  }, error = function(e) {
    message(paste("  HMM fitting failed for", n_states, "states:", e$message))
    return(NULL)
  })
}

# In-sample performance for REM
get_insample_perf_MLE <- function(coeff, stats, reh) {
  # Compute log-likelihood based predictions
  n_events <- reh$M
  n_actors <- reh$N
  
  # Get risk set
  riskset <- reh$riskset
  
  # Placeholder for actual implementation - return basic metrics
  list(
    n_events = n_events,
    n_actors = n_actors
  )
}

# Out-of-sample performance for REM
get_outofsample_perf_MLE <- function(coeff, stats, reh, M_train, M_test) {
  # Placeholder - compute predictive performance metrics
  list(
    M_train = M_train,
    M_test = M_test,
    coeff_count = length(coeff)
  )
}

# -----------------------------------------------------------------------------
# 4. Load Apollo 13 Data
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 70), "\n")
cat("  APOLLO 13 HMMREM ANALYSIS\n")
cat("  With Comprehensive Metrics\n")
cat(strrep("=", 70), "\n\n")

cat("Loading Apollo 13 data...\n")
load("Apollo13_completeData.RData")
time_sender_receiver <- Merged_ALL_parts_Apollo

cat("Data loaded:", nrow(time_sender_receiver), "events\n")
cat("Actors:", length(unique(c(time_sender_receiver$sender, time_sender_receiver$receiver))), "\n")

# -----------------------------------------------------------------------------
# 5. Compute Adjusted Time Differences
# -----------------------------------------------------------------------------

cat("\nComputing adjusted time differences...\n")

time_differences <- diff(time_sender_receiver$time)
time_sender_receiver$adjusted_differences <- c(0, time_differences)

time_sender_receiver <- time_sender_receiver %>%
  group_by(time) %>%
  mutate(adjusted_differences = first(adjusted_differences)) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 6. HMM Model Comparison (1, 2, 3, 4 States)
# -----------------------------------------------------------------------------

cat("\n", strrep("-", 70), "\n")
cat("  HMM MODEL COMPARISON\n")
cat(strrep("-", 70), "\n\n")

# Prepare edgelist
edgelist <- time_sender_receiver
edgelist$message <- NULL
edgelist <- na.omit(edgelist)
edgelist$adjusted_differences[1] <- 1

# Prepare data for HMM
dt <- edgelist$adjusted_differences
HMMdf <- data.frame(
  ID = rep("a", length(dt)),
  step = dt,
  angle = NA
)

hmm_data <- prepData(data = HMMdf, coordNames = NULL)

# Store results for all state counts
hmm_results <- list()

for (n_states in states_to_fit) {
  cat("Fitting HMM with", n_states, "state(s)...\n")
  
  hmm_fit <- fit_hmm_safe(hmm_data, n_states)
  
  if (!is.null(hmm_fit)) {
    bic <- AIC(hmm_fit, k = log(nrow(hmm_data)))
    aic <- AIC(hmm_fit)
    
    if (n_states == 1) {
      conf_metrics <- list(avg_confidence = 1.0, pct_uncertain_tau1 = 0, 
                           pct_uncertain_tau2 = 0, state_probs = NULL, max_probs = NULL)
      viterbi_states <- rep(1, nrow(edgelist))
    } else {
      conf_metrics <- compute_confidence_metrics(hmm_fit, tau1, tau2)
      viterbi_states <- momentuHMM::viterbi(hmm_fit)
    }
    
    hmm_results[[paste0("m", n_states)]] <- list(
      n_states = n_states,
      bic = bic,
      aic = aic,
      avg_confidence = conf_metrics$avg_confidence,
      pct_uncertain_tau1 = conf_metrics$pct_uncertain_tau1,
      pct_uncertain_tau2 = conf_metrics$pct_uncertain_tau2,
      state_probs = conf_metrics$state_probs,
      max_probs = conf_metrics$max_probs,
      viterbi_states = viterbi_states,
      step_params = hmm_fit$mle$step,
      gamma = if (n_states > 1) hmm_fit$mle$gamma else NULL,
      fit = hmm_fit
    )
    
    cat("  BIC:", round(bic, 2), "\n")
    cat("  AIC:", round(aic, 2), "\n")
    cat("  Avg Confidence:", round(conf_metrics$avg_confidence * 100, 2), "%\n")
    cat("  % Uncertain (tau <", tau1, "):", round(conf_metrics$pct_uncertain_tau1, 2), "%\n")
    cat("  % Uncertain (tau <", tau2, "):", round(conf_metrics$pct_uncertain_tau2, 2), "%\n")
    cat("\n")
  }
}

# Create summary table
hmm_summary <- data.frame(
  States = sapply(hmm_results, function(x) x$n_states),
  BIC = sapply(hmm_results, function(x) x$bic),
  AIC = sapply(hmm_results, function(x) x$aic),
  Confidence = sapply(hmm_results, function(x) x$avg_confidence * 100),
  Uncertain_tau1 = sapply(hmm_results, function(x) x$pct_uncertain_tau1),
  Uncertain_tau2 = sapply(hmm_results, function(x) x$pct_uncertain_tau2)
)

cat("\n--- HMM Model Selection Summary ---\n\n")
print(knitr::kable(hmm_summary, digits = 2, 
                   col.names = c("States", "BIC", "AIC", 
                                 "Confidence (%)", 
                                 paste0("Uncertain (tau<", tau1, ")"), 
                                 paste0("Uncertain (tau<", tau2, ")"))))

# Best model by BIC
best_model <- names(hmm_results)[which.min(sapply(hmm_results, function(x) x$bic))]
cat("\nOptimal number of states (by BIC):", hmm_results[[best_model]]$n_states, "\n")

# Display transition matrices for all models
cat("\n", strrep("-", 70), "\n")
cat("  TRANSITION PROBABILITY MATRICES\n")
cat(strrep("-", 70), "\n\n")

for (model_name in names(hmm_results)) {
  result <- hmm_results[[model_name]]
  if (result$n_states > 1) {
    cat("--- Transition Matrix (", result$n_states, "-State Model) ---\n", sep = "")
    trans_mat <- result$fit$mle$gamma
    rownames(trans_mat) <- paste0("From State ", 1:result$n_states)
    colnames(trans_mat) <- paste0("To State ", 1:result$n_states)
    print(knitr::kable(trans_mat, digits = 4))
    cat("\n")
  }
}

# -----------------------------------------------------------------------------
# 7. Use Best HMM Model for State Prediction
# -----------------------------------------------------------------------------

# Use 2-state model (as in original analysis)
nstates <- 2
best_fit <- hmm_results[["m2"]]$fit
state_sequence <- hmm_results[["m2"]]$viterbi_states
state_probs <- hmm_results[["m2"]]$state_probs

edgelist$Predicted <- state_sequence
state_counts <- table(edgelist$Predicted)

cat("\n--- State Distribution (2-State Model) ---\n")
print(state_counts)

# Prepare segmentation
edgelist$row_index <- seq_len(nrow(edgelist))
means_by_group <- aggregate(adjusted_differences ~ Predicted, data = edgelist, FUN = mean)
sorted_index <- order(means_by_group$adjusted_differences, decreasing = TRUE)

edgelist <- edgelist %>%
  mutate(row_index = row_number(),
         segment_id = cumsum(c(TRUE, diff(Predicted) != 0)))

rects <- edgelist %>%
  group_by(segment_id, Predicted) %>%
  summarize(start = min(row_index), end = max(row_index) + 1, .groups = 'drop') %>%
  mutate(Color = case_when(
    Predicted == as.numeric(sorted_index[1]) ~ "Low",
    Predicted == as.numeric(sorted_index[2]) ~ "High"
  ))

edgelist$column0 <- ifelse(edgelist$Predicted == sorted_index[2], 1, 0)

# -----------------------------------------------------------------------------
# 8. Confidence Analysis and Visualization
# -----------------------------------------------------------------------------

cat("\n", strrep("-", 70), "\n")
cat("  CONFIDENCE ANALYSIS (2-State Model)\n")
cat(strrep("-", 70), "\n\n")

# Add confidence scores to edgelist
if (!is.null(state_probs)) {
  edgelist$state_prob <- apply(state_probs, 1, max)
  edgelist$confident_tau1 <- edgelist$state_prob >= tau1
  edgelist$confident_tau2 <- edgelist$state_prob >= tau2
  
  cat("Confidence Statistics:\n")
  cat("  Mean confidence:", round(mean(edgelist$state_prob) * 100, 2), "%\n")
  cat("  Median confidence:", round(median(edgelist$state_prob) * 100, 2), "%\n")
  cat("  Min confidence:", round(min(edgelist$state_prob) * 100, 2), "%\n")
  cat("  Max confidence:", round(max(edgelist$state_prob) * 100, 2), "%\n")
  cat("\n")
  cat("  Events with confidence >= ", tau1, ": ", sum(edgelist$confident_tau1), 
      " (", round(mean(edgelist$confident_tau1) * 100, 1), "%)\n", sep = "")
  cat("  Events with confidence >= ", tau2, ": ", sum(edgelist$confident_tau2), 
      " (", round(mean(edgelist$confident_tau2) * 100, 1), "%)\n", sep = "")
  
  # Confidence distribution by state
  cat("\nConfidence by State:\n")
  conf_by_state <- edgelist %>%
    group_by(Predicted) %>%
    summarize(
      n = n(),
      mean_conf = mean(state_prob),
      sd_conf = sd(state_prob),
      pct_confident_tau1 = mean(confident_tau1) * 100,
      pct_confident_tau2 = mean(confident_tau2) * 100
    )
  print(knitr::kable(conf_by_state, digits = 3))
}

# Plot: Confidence over time
p_confidence <- ggplot(edgelist[1:1000, ], aes(x = row_index, y = state_prob * 100)) +
  geom_line(color = "#3498DB", alpha = 0.7) +
  geom_hline(yintercept = tau1 * 100, color = "#E74C3C", linetype = "dashed", size = 0.8) +
  geom_hline(yintercept = tau2 * 100, color = "#F39C12", linetype = "dashed", size = 0.8) +
  labs(title = "State Classification Confidence Over Time",
       subtitle = paste0("Red dashed: tau1=", tau1, ", Orange dashed: tau2=", tau2),
       x = "Event Index", y = "Confidence (%)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

print(p_confidence)
ggsave("Apollo_confidence_over_time.png", p_confidence, width = 10, height = 5, dpi = 300)
cat("\nSaved: Apollo_confidence_over_time.png\n")

# Plot: HMM segmentation
PTimeState <- ggplot() + 
  geom_rect(data = rects[1:22, ], 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = Color), alpha = 1) + 
  scale_fill_manual(values = c("Low" = "#5DADE2", "High" = "#FFC0CB"), 
                    limits = c("Low", "High")) +
  geom_line(data = edgelist[1:1000, ], 
            aes(x = row_index, y = adjusted_differences), size = 0.5) +  
  geom_vline(xintercept = 249, color = "black", linetype = "dashed", size = 1) + 
  labs(x = "Row Index", y = "Event Frequency (Δt)", 
       title = "HMM State Segmentation (Apollo 13)") + 
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(colour = "#00000080", fill = NA, size = 1.5),
    text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank()
  )

print(PTimeState)
ggsave("Apollo_HMM_segmentation.png", PTimeState, width = 12, height = 6, dpi = 300)
cat("Saved: Apollo_HMM_segmentation.png\n")

# -----------------------------------------------------------------------------
# 9. Prepare Covariates for REM
# -----------------------------------------------------------------------------

cat("\n", strrep("-", 70), "\n")
cat("  RELATIONAL EVENT MODELS\n")
cat(strrep("-", 70), "\n\n")

names(edgelist)[2:3] <- c("actor1", "actor2")
actors <- sort(unique(c(edgelist$actor1, edgelist$actor2)))

info <- data.frame(
  id = actors,
  time = 0,
  team = ifelse(actors %in% c("CDR", "LMP", "CMP"), "air", "ground")
)

info <- dummy_cols(info, select_columns = "id")
info <- dummy_cols(info, select_columns = "team")

getTie <- function(var1, var2, actors) {
  x <- var1 %*% t(var2)
  rownames(x) <- colnames(x) <- actors
  x
}

air_to_air <- getTie(info$team_air, info$team_air, info$id)
ground_to_CAPCOM <- getTie(info$team_ground, info$id_CAPCOM, info$id)
air_to_CAPCOM <- getTie(info$team_air, info$id_CAPCOM, info$id)
CAPCOM_to_air <- getTie(info$id_CAPCOM, info$team_air, info$id)
ground_to_FLIGHT <- getTie(info$team_ground, info$id_FLIGHT, info$id)
FLIGHT_to_ground <- getTie(info$id_FLIGHT, info$team_ground, info$id)

# Store REM results
rem_results <- list()

# -----------------------------------------------------------------------------
# 10. Model 1: Simple REM (No State Information)
# -----------------------------------------------------------------------------

cat("Fitting Model 1: Simple REM...\n")

stats <- ~ 1 + 
  rrankSend() + 
  psABBA() + psABBY() + psABXA() + psABAY() +
  reciprocity(scaling = "prop") +
  outdegreeSender(scaling = "prop") +
  tie(variable = "air_to_air", attr_dyads = air_to_air) +
  tie(variable = "ground_to_CAPCOM", attr_dyads = ground_to_CAPCOM) +
  tie(variable = "ground_to_FLIGHT", attr_dyads = ground_to_FLIGHT) +
  tie(variable = "FLIGHT_to_ground", attr_dyads = FLIGHT_to_ground) +
  tie(variable = "air_to_CAPCOM", attr_dyads = air_to_CAPCOM) +
  (tie(variable = "CAPCOM_to_air", attr_dyads = CAPCOM_to_air)) :
  (rrankReceive())

# Split training and test sets
edgelist_train <- edgelist[1:N_train, ]
edgelist_test  <- edgelist[(N_train + 1):(N_train + N_test), ]

# Create remify object for training data
reh_train <- remify::remify(
  edgelist = edgelist_train,
  riskset = "active",
  model = "tie",
  directed = TRUE,
  origin = 0
)

# Compute statistics
out_train <- remstats(
  reh = reh_train,
  tie_effects = stats,
  attr_actors = info
)

# Fit REM
fit_train <- remstimate::remstimate(
  reh = reh_train,
  stats = out_train,
  method = "MLE"
)

summary1 <- summary(fit_train)
cat("\n--- Model 1: Simple REM ---\n")
print(summary1)

rem_results[["simple"]] <- list(
  model = "Simple REM",
  coefs = summary1$coefsTab,
  bic = summary1$BIC,
  aic = summary1$AIC
)

# -----------------------------------------------------------------------------
# 11. Model 2: REM with State Effect (HMMREM)
# -----------------------------------------------------------------------------

cat("\nFitting Model 2: REM + State Effect (HMMREM)...\n")

stats2train <- ~ 1 + 
  rrankSend() + 
  psABBA() + psABBY() + psABXA() + psABAY() +
  reciprocity(scaling = "prop") +
  outdegreeSender(scaling = "prop") +
  tie(variable = "air_to_air", attr_dyads = air_to_air) +
  tie(variable = "ground_to_CAPCOM", attr_dyads = ground_to_CAPCOM) +
  tie(variable = "ground_to_FLIGHT", attr_dyads = ground_to_FLIGHT) +
  tie(variable = "FLIGHT_to_ground", attr_dyads = FLIGHT_to_ground) +
  tie(variable = "air_to_CAPCOM", attr_dyads = air_to_CAPCOM) +
  (tie(variable = "CAPCOM_to_air", attr_dyads = CAPCOM_to_air)) :
  (rrankReceive() + event(x = edgelist_train$column0, "PredictedState1")) +
  rrankReceive() : event(x = edgelist_train$column0, "PredictedState1")

out_train2 <- remstats(
  reh = reh_train,
  tie_effects = stats2train,
  attr_actors = info
)

fit_train2 <- remstimate::remstimate(
  reh = reh_train,
  stats = out_train2,
  method = "MLE"
)

summary2 <- summary(fit_train2)
cat("\n--- Model 2: REM + State Effect ---\n")
print(summary2)

rem_results[["state"]] <- list(
  model = "REM + State",
  coefs = summary2$coefsTab,
  bic = summary2$BIC,
  aic = summary2$AIC
)

# -----------------------------------------------------------------------------
# 12. Model Comparison
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 70), "\n")
cat("  MODEL COMPARISON SUMMARY\n")
cat(strrep("=", 70), "\n\n")

# HMM comparison
cat("--- HMM State Selection ---\n")
print(knitr::kable(hmm_summary, digits = 2))

# REM comparison
cat("\n--- REM Model Comparison ---\n\n")

rem_comparison <- data.frame(
  Model = sapply(rem_results, function(x) x$model),
  BIC = sapply(rem_results, function(x) x$bic),
  AIC = sapply(rem_results, function(x) x$aic)
)

print(knitr::kable(rem_comparison, digits = 2))

# Delta BIC
rem_comparison$Delta_BIC <- rem_comparison$BIC - min(rem_comparison$BIC)

cat("\n--- Delta BIC (relative to best model) ---\n")
print(knitr::kable(rem_comparison[, c("Model", "BIC", "Delta_BIC")], digits = 2))

# -----------------------------------------------------------------------------
# 13. Final Summary Table
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 70), "\n")
cat("  FINAL RESULTS\n")
cat(strrep("=", 70), "\n\n")

cat("Best HMM Model:\n")
cat("  Number of states:", hmm_results[[best_model]]$n_states, "\n")
cat("  BIC:", round(hmm_results[[best_model]]$bic, 2), "\n")
cat("  Avg Confidence:", round(hmm_results[[best_model]]$avg_confidence * 100, 2), "%\n")
cat("  % Uncertain (tau <", tau1, "):", round(hmm_results[[best_model]]$pct_uncertain_tau1, 2), "%\n")
cat("  % Uncertain (tau <", tau2, "):", round(hmm_results[[best_model]]$pct_uncertain_tau2, 2), "%\n")

best_rem <- names(rem_results)[which.min(sapply(rem_results, function(x) x$bic))]
cat("\nBest REM Model:\n")
cat("  Model:", rem_results[[best_rem]]$model, "\n")
cat("  BIC:", round(rem_results[[best_rem]]$bic, 2), "\n")
cat("  AIC:", round(rem_results[[best_rem]]$aic, 2), "\n")

# Save results for export
assign("hmm_results_apollo", hmm_results, envir = .GlobalEnv)
assign("rem_results_apollo", rem_results, envir = .GlobalEnv)
assign("hmm_summary_apollo", hmm_summary, envir = .GlobalEnv)
assign("edgelist_apollo", edgelist, envir = .GlobalEnv)

# -----------------------------------------------------------------------------
# 14. Export Results to Excel
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 70), "\n")
cat("  EXPORTING RESULTS TO EXCEL\n")
cat(strrep("=", 70), "\n\n")

# Check if openxlsx is available, if not install it
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  message("Installing package: openxlsx")
  install.packages("openxlsx", repos = "https://cloud.r-project.org")
}
library(openxlsx)

# Create workbook
wb <- createWorkbook()

# Sheet 1: HMM Model Comparison
addWorksheet(wb, "HMM_Model_Comparison")
writeData(wb, "HMM_Model_Comparison", hmm_summary, startRow = 1, startCol = 1)

# Sheet 2: Transition Matrices
addWorksheet(wb, "Transition_Matrices")
row_offset <- 1
for (model_name in names(hmm_results)) {
  result <- hmm_results[[model_name]]
  if (result$n_states > 1) {
    # Write title
    writeData(wb, "Transition_Matrices", 
              paste0(result$n_states, "-State Transition Matrix"),
              startRow = row_offset, startCol = 1)
    row_offset <- row_offset + 1
    
    # Write matrix
    trans_mat <- result$fit$mle$gamma
    trans_df <- as.data.frame(trans_mat)
    rownames(trans_df) <- paste0("From State ", 1:result$n_states)
    colnames(trans_df) <- paste0("To State ", 1:result$n_states)
    trans_df <- cbind(State = rownames(trans_df), trans_df)
    writeData(wb, "Transition_Matrices", trans_df, 
              startRow = row_offset, startCol = 1)
    row_offset <- row_offset + nrow(trans_df) + 2
  }
}

# Sheet 3: REM Model Comparison
addWorksheet(wb, "REM_Model_Comparison")
writeData(wb, "REM_Model_Comparison", rem_comparison, startRow = 1, startCol = 1)

# Sheet 4: Simple REM Coefficients
if (!is.null(rem_results$simple$coefs)) {
  addWorksheet(wb, "Simple_REM_Coefficients")
  coefs_simple <- as.data.frame(rem_results$simple$coefs)
  coefs_simple <- cbind(Parameter = rownames(coefs_simple), coefs_simple)
  writeData(wb, "Simple_REM_Coefficients", coefs_simple, startRow = 1, startCol = 1)
}

# Sheet 5: HMMREM Coefficients
if (!is.null(rem_results$state$coefs)) {
  addWorksheet(wb, "HMMREM_Coefficients")
  coefs_state <- as.data.frame(rem_results$state$coefs)
  coefs_state <- cbind(Parameter = rownames(coefs_state), coefs_state)
  writeData(wb, "HMMREM_Coefficients", coefs_state, startRow = 1, startCol = 1)
}

# Sheet 6: Confidence Statistics (2-State Model)
addWorksheet(wb, "Confidence_Statistics")

# Extract state probabilities correctly
state_probs_max <- as.numeric(hmm_results$m2$state_probs_max)
n_events_total <- length(state_probs_max)
n_conf_tau1 <- sum(state_probs_max >= tau1, na.rm = TRUE)
n_conf_tau2 <- sum(state_probs_max >= tau2, na.rm = TRUE)

conf_stats <- data.frame(
  Metric = c(
    "Mean Confidence", 
    "Median Confidence", 
    "Min Confidence", 
    "Max Confidence",
    paste0("Count: Confidence >= ", tau1),
    paste0("Count: Confidence >= ", tau2)
  ),
  Value = c(
    mean(state_probs_max, na.rm = TRUE),
    median(state_probs_max, na.rm = TRUE),
    min(state_probs_max, na.rm = TRUE),
    max(state_probs_max, na.rm = TRUE),
    n_conf_tau1,
    n_conf_tau2
  ),
  Percentage = c(
    mean(state_probs_max, na.rm = TRUE) * 100,
    median(state_probs_max, na.rm = TRUE) * 100,
    min(state_probs_max, na.rm = TRUE) * 100,
    max(state_probs_max, na.rm = TRUE) * 100,
    (n_conf_tau1 / n_events_total) * 100,
    (n_conf_tau2 / n_events_total) * 100
  )
)
writeData(wb, "Confidence_Statistics", conf_stats, startRow = 1, startCol = 1)

# Sheet 7: Confidence by State
addWorksheet(wb, "Confidence_by_State")

# Create a simpler confidence by state summary
conf_by_state_list <- list()
viterbi_states <- hmm_results$m2$viterbi_states
for (state_id in sort(unique(viterbi_states))) {
  state_indices <- which(viterbi_states == state_id)
  state_conf <- state_probs_max[state_indices]
  
  conf_by_state_list[[as.character(state_id)]] <- data.frame(
    Predicted_State = state_id,
    n = length(state_indices),
    mean_conf = mean(state_conf, na.rm = TRUE),
    sd_conf = sd(state_conf, na.rm = TRUE),
    pct_confident_tau1 = mean(state_conf >= tau1, na.rm = TRUE) * 100,
    pct_confident_tau2 = mean(state_conf >= tau2, na.rm = TRUE) * 100
  )
}
conf_by_state <- do.call(rbind, conf_by_state_list)
rownames(conf_by_state) <- NULL
writeData(wb, "Confidence_by_State", conf_by_state, startRow = 1, startCol = 1)

# Save workbook
excel_filename <- "Apollo_HMMREM_Results.xlsx"
saveWorkbook(wb, excel_filename, overwrite = TRUE)

cat("Excel file created:", excel_filename, "\n")
cat("Sheets included:\n")
cat("  1. HMM_Model_Comparison\n")
cat("  2. Transition_Matrices\n")
cat("  3. REM_Model_Comparison\n")
cat("  4. Simple_REM_Coefficients\n")
cat("  5. HMMREM_Coefficients\n")
cat("  6. Confidence_Statistics\n")
cat("  7. Confidence_by_State\n")

cat("\n\nAnalysis complete!\n")
cat("Results saved to global environment: hmm_results_apollo, rem_results_apollo, hmm_summary_apollo, edgelist_apollo\n")
cat("Excel report saved to:", excel_filename, "\n")
