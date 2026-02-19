# =============================================================================
# HMMREM Scenario Analysis: EASY, MEDIUM, HARD
# Comparing REM with 1, 2, 3, 4 hidden states across difficulty levels
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Setup and Package Loading
# -----------------------------------------------------------------------------

required_packages <- c("remstats", "remstimate", "dplyr", "ggplot2", 
                       "plotly", "momentuHMM", "remulate", "tidyr", "knitr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed."))
  } else {
    library(pkg, character.only = TRUE)
  }
}

set.seed(1234)

# -----------------------------------------------------------------------------
# 2. Define Scenario Parameters
# -----------------------------------------------------------------------------

scenarios <- list(
  EASY = list(
    name = "EASY",
    n_events = 3000,
    transition_probs = matrix(c(0.99, 0.01,
                                0.02, 0.98), nrow = 2, byrow = TRUE),
    baseline1 = -8,
    baseline2 = -6.5,
    outdegree1 = 0.2,
    outdegree2 = 0.5,
    inertia1 = 0.2,
    inertia2 = 0.4
  ),
  MEDIUM = list(
    name = "MEDIUM",
    n_events = 3000,
    transition_probs = matrix(c(0.95, 0.05,
                                0.10, 0.90), nrow = 2, byrow = TRUE),
    baseline1 = -8,
    baseline2 = -6.5,
    outdegree1 = 0.2,
    outdegree2 = 0.5,
    inertia1 = 0.2,
    inertia2 = 0.4
  ),
  HARD = list(
    name = "HARD",
    n_events = 3000,
    transition_probs = matrix(c(0.95, 0.05,
                                0.10, 0.90), nrow = 2, byrow = TRUE),
    baseline1 = -6,
    baseline2 = -5,
    outdegree1 = 0.2,
    outdegree2 = 0.45,
    inertia1 = 0.25,
    inertia2 = 0.4
  ),
    EHARD = list(
    name = "EHARD",
    n_events = 3000,
    transition_probs = matrix(c(0.97, 0.03,
                                0.26, 0.74), nrow = 2, byrow = TRUE),
    baseline1 = -6,
    baseline2 = -5,
    outdegree1 = 0.2,
    outdegree2 = 0.45,
    inertia1 = 0.25,
    inertia2 = 0.4
  )
)

# Number of replications
R <- 3

# Number of actors
n_actors <- 15

# States to test in HMM fitting
states_to_fit <- 1:4

# Confidence thresholds
tau1 <- 0.9
tau2 <- 0.8

# -----------------------------------------------------------------------------
# 3. Helper Functions
# -----------------------------------------------------------------------------

simulate_hidden_states <- function(n_events, transition_probs, initial_probs = c(1, 0)) {
  m <- nrow(transition_probs)
  hidden_states <- numeric(n_events)
  
  for (t in seq_len(n_events)) {
    if (t == 1) {
      hidden_states[t] <- sample(m, 1, prob = initial_probs)
    } else {
      hidden_states[t] <- sample(m, 1, prob = transition_probs[hidden_states[t - 1], ])
    }
  }
  return(hidden_states)
}

check_valid_states <- function(hidden_states, min_segment_length = 1) {
  # Check if all state segments have at least min_segment_length
  # For harder scenarios, we relax this constraint
  state_runs <- rle(hidden_states)
  min_len <- min(state_runs$lengths)
  return(min_len >= min_segment_length)
}

generate_valid_hidden_states <- function(n_events, transition_probs, initial_probs = c(1, 0), 
                                          max_attempts = 100, min_segment_length = 2) {
  # Try to generate hidden states with minimum segment length

  # If we can't after max_attempts, relax the constraint
  
  for (attempt in 1:max_attempts) {
    hidden_states <- simulate_hidden_states(n_events, transition_probs, initial_probs)
    
    if (check_valid_states(hidden_states, min_segment_length)) {
      return(hidden_states)
    }
    
    # Every 20 attempts, print a progress message
    if (attempt %% 20 == 0) {
      cat("    Attempt", attempt, "to generate valid states...\n")
    }
  }
  
  # If we couldn't get valid states with strict constraint, relax it
  cat("    Could not generate states with min segment length", min_segment_length, 
      "- relaxing constraint\n")
  
  # Just return the last generated states - they may have some length-1 segments
  # but that's acceptable for harder scenarios
  return(hidden_states)
}

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
      pct_uncertain_tau2 = pct_uncertain_tau2
    ))
  }, error = function(e) {
    return(list(avg_confidence = NA, pct_uncertain_tau1 = NA, pct_uncertain_tau2 = NA))
  })
}

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

calculate_accuracy <- function(true_states, predicted_states, n_states) {
  if (n_states != 2) return(NA)
  acc1 <- mean(true_states == predicted_states) * 100
  predicted_flipped <- ifelse(predicted_states == 1, 2, 1)
  acc2 <- mean(true_states == predicted_flipped) * 100
  return(max(acc1, acc2))
}

# -----------------------------------------------------------------------------
# 4. Aggregation Functions
# -----------------------------------------------------------------------------

aggregate_hmm_results <- function(scenario_name, replications) {
  n_rep <- length(replications)
  state_counts <- names(replications[[1]]$hmm_results)
  agg_results <- data.frame()
  
  for (sc in state_counts) {
    m <- as.numeric(gsub("m", "", sc))
    bics <- sapply(1:n_rep, function(i) replications[[i]]$hmm_results[[sc]]$bic)
    aics <- sapply(1:n_rep, function(i) replications[[i]]$hmm_results[[sc]]$aic)
    confidences <- sapply(1:n_rep, function(i) replications[[i]]$hmm_results[[sc]]$avg_confidence)
    uncertainties_tau1 <- sapply(1:n_rep, function(i) replications[[i]]$hmm_results[[sc]]$pct_uncertain_tau1)
    uncertainties_tau2 <- sapply(1:n_rep, function(i) replications[[i]]$hmm_results[[sc]]$pct_uncertain_tau2)
    accuracies <- sapply(1:n_rep, function(i) replications[[i]]$hmm_results[[sc]]$accuracy)
    
    row <- data.frame(
      Scenario = scenario_name, States = m,
      BIC_mean = mean(bics, na.rm = TRUE), BIC_sd = sd(bics, na.rm = TRUE),
      AIC_mean = mean(aics, na.rm = TRUE), AIC_sd = sd(aics, na.rm = TRUE),
      Confidence_mean = mean(confidences, na.rm = TRUE),
      Confidence_sd = sd(confidences, na.rm = TRUE),
      Uncertainty_tau1_mean = mean(uncertainties_tau1, na.rm = TRUE),
      Uncertainty_tau1_sd = sd(uncertainties_tau1, na.rm = TRUE),
      Uncertainty_tau2_mean = mean(uncertainties_tau2, na.rm = TRUE),
      Uncertainty_tau2_sd = sd(uncertainties_tau2, na.rm = TRUE),
      Accuracy_mean = mean(accuracies, na.rm = TRUE),
      Accuracy_sd = sd(accuracies, na.rm = TRUE)
    )
    agg_results <- rbind(agg_results, row)
  }
  return(agg_results)
}

aggregate_rem_results <- function(replications) {
  n_rep <- length(replications)
  if (length(replications[[1]]$rem_results) == 0) return(NULL)
  
  model_types <- names(replications[[1]]$rem_results)
  all_rem_results <- list()
  
  for (model in model_types) {
    all_coefs <- lapply(1:n_rep, function(i) {
      coef_tab <- replications[[i]]$rem_results[[model]]$coefs
      if (!is.null(coef_tab)) {
        # Handle different column name formats
        est_col <- if ("Estimate" %in% colnames(coef_tab)) "Estimate" else 1
        se_col <- if ("Std. Error" %in% colnames(coef_tab)) "Std. Error" else 
                  if ("Std.Err" %in% colnames(coef_tab)) "Std.Err" else 2
        z_col <- if ("z value" %in% colnames(coef_tab)) "z value" else 3
        
        # Get z-values and compute p-values properly (two-tailed test)
        z_values <- as.numeric(coef_tab[, z_col])
        # Compute p-values from z-values using normal distribution
        pvalues <- 2 * pnorm(-abs(z_values))
        # Set minimum p-value to avoid log(0)
        pvalues[pvalues < 1e-300] <- 1e-300
        
        data.frame(Variable = rownames(coef_tab), 
                   Estimate = as.numeric(coef_tab[, est_col]),
                   SE = as.numeric(coef_tab[, se_col]),
                   z_value = z_values,
                   pvalue = pvalues)
      }
    })
    
    combined <- do.call(rbind, all_coefs)
    
    agg <- combined %>%
      group_by(Variable) %>%
      summarize(
        Estimate_mean = mean(Estimate, na.rm = TRUE),
        Estimate_sd = sd(Estimate, na.rm = TRUE),
        SE_mean = mean(SE, na.rm = TRUE),
        z_value_mean = mean(z_value, na.rm = TRUE),
        z_value_sd = sd(z_value, na.rm = TRUE),
        pvalue_mean = mean(pvalue, na.rm = TRUE),
        pvalue_sd = sd(pvalue, na.rm = TRUE),
        log10_pvalue_mean = mean(-log10(pvalue), na.rm = TRUE),
        sig_rate = mean(pvalue < 0.05, na.rm = TRUE) * 100,
        .groups = 'drop'
      )
    
    bics <- sapply(1:n_rep, function(i) replications[[i]]$rem_results[[model]]$bic)
    
    all_rem_results[[model]] <- list(
      model = replications[[1]]$rem_results[[model]]$model,
      coefficients = agg,
      bic_mean = mean(bics, na.rm = TRUE),
      bic_sd = sd(bics, na.rm = TRUE)
    )
  }
  return(all_rem_results)
}

# -----------------------------------------------------------------------------
# 5. Visualization Functions
# -----------------------------------------------------------------------------

create_bic_comparison_plot <- function(all_hmm_results, R_val) {
  combined <- do.call(rbind, all_hmm_results)
  combined$Scenario <- factor(combined$Scenario, levels = c("EASY", "MEDIUM", "HARD", "EHARD"))
  combined$States <- factor(combined$States)
  
  p <- ggplot(combined, aes(x = States, y = BIC_mean, fill = Scenario)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = BIC_mean - BIC_sd, ymax = BIC_mean + BIC_sd),
                  position = position_dodge(width = 0.8), width = 0.2) +
    scale_fill_manual(values = c("EASY" = "#2ECC71", "MEDIUM" = "#F1C40F", "HARD" = "#E74C3C", "EHARD" = "#9B59B6")) +
    labs(title = "HMM Model Comparison: BIC by Number of States",
         subtitle = paste0("Based on ", R_val, " replications per scenario"),
         x = "Number of Hidden States", y = "BIC (mean +/- SD)",
         fill = "Scenario Difficulty") +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.position = "bottom")
  
  return(p)
}

create_confidence_plot <- function(all_hmm_results, tau1_val, tau2_val) {
  combined <- do.call(rbind, all_hmm_results)
  combined$Scenario <- factor(combined$Scenario, levels = c("EASY", "MEDIUM", "HARD", "EHARD"))
  combined_2state <- combined %>% filter(States == 2)
  
  # Create metric labels for consistent naming
  metric_tau1 <- paste0("% Uncertain (tau=", tau1_val, ")")
  metric_tau2 <- paste0("% Uncertain (tau=", tau2_val, ")")
  
  plot_data <- combined_2state %>%
    select(Scenario, Confidence_mean, Uncertainty_tau1_mean, Uncertainty_tau2_mean) %>%
    mutate(Confidence_pct = Confidence_mean * 100) %>%
    pivot_longer(cols = c(Confidence_pct, Uncertainty_tau1_mean, Uncertainty_tau2_mean),
                 names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = case_when(
      Metric == "Confidence_pct" ~ "Avg Confidence",
      Metric == "Uncertainty_tau1_mean" ~ metric_tau1,
      Metric == "Uncertainty_tau2_mean" ~ metric_tau2
    ))
  
  # Create color mapping
  color_map <- c("Avg Confidence" = "#3498DB", 
                 metric_tau1 = "#E74C3C",
                 metric_tau2 = "#F39C12")
  
  p <- ggplot(plot_data, aes(x = Scenario, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = color_map) +
    labs(title = "State Classification Confidence (2-State HMM)",
         subtitle = paste0("Uncertainty thresholds: tau1=", tau1_val, ", tau2=", tau2_val),
         x = "Scenario Difficulty", y = "Percentage (%)", fill = "Metric") +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold"), legend.position = "bottom")
  
  return(p)
}

# -----------------------------------------------------------------------------
# 6. Main Execution (Inline like the original working code)
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 70), "\n")
cat("  HMMREM SCENARIO ANALYSIS\n")
cat("  Difficulty Levels: EASY, MEDIUM, HARD, EHARD\n")
cat("  Replications:", R, "\n")
cat("  States to fit: 1, 2, 3, 4\n")
cat(strrep("=", 70), "\n")

all_results <- list()
all_hmm_agg <- list()
all_rem_agg <- list()

# Loop through scenarios
for (scenario_name in names(scenarios)) {
  scenario <- scenarios[[scenario_name]]
  
  cat("\n", strrep("=", 60), "\n")
  cat("Running Scenario:", scenario_name, "\n")
  cat("n_events =", scenario$n_events, "\n")
  cat(strrep("=", 60), "\n")
  
  replications <- list()
  
  for (iter in 1:R) {
    cat("  Replication", iter, "\n")
    
    # Create actor attributes - assign to global env for formula capture
    actors <- 1:n_actors
    attr_actors <<- data.frame(
      name = actors,
      time = rep(0, n_actors),
      sex = sample(0:1, n_actors, replace = TRUE),
      age = sample(0:1, n_actors, replace = TRUE)
    )
    
    # Simulate hidden states with timeout protection
    # For EASY scenario, we require min segment length of 2
    # For MEDIUM/HARD, we relax this to allow length-1 segments if needed
    min_seg_len <- if (scenario_name == "EASY") 2 else 1
    hidden_states <- generate_valid_hidden_states(
      scenario$n_events, 
      scenario$transition_probs,
      max_attempts = 50,
      min_segment_length = min_seg_len
    )
    
    state_runs <- rle(hidden_states)
    state_lengths <- state_runs$lengths
    state_values <- state_runs$values
    
    # Build effects based on scenario - using global attr_actors
    if (scenario_name == "EASY") {
      effects1 <- ~
        remulate::baseline(-8) +
        remulate::difference(0.3, "sex", attr_actors, scaling = "std") +
        remulate::difference(0.4, "age", attr_actors, scaling = "std") +
        remulate::outdegreeReceiver(0.2, scaling = "std") +
        remulate::inertia(0.2, scaling = "std")
      
      effects2 <- ~
        remulate::baseline(-6.5) +
        remulate::difference(0.3, "sex", attr_actors, scaling = "std") +
        remulate::difference(0.4, "age", attr_actors, scaling = "std") +
        remulate::outdegreeReceiver(0.5, scaling = "std") +
        remulate::inertia(0.4, scaling = "std")
        
    } else if (scenario_name == "MEDIUM") {
      effects1 <- ~
        remulate::baseline(-8) +
        remulate::difference(0.3, "sex", attr_actors, scaling = "std") +
        remulate::difference(0.4, "age", attr_actors, scaling = "std") +
        remulate::outdegreeReceiver(0.2, scaling = "std") +
        remulate::inertia(0.2, scaling = "std")
      
      effects2 <- ~
        remulate::baseline(-6.5) +
        remulate::difference(0.3, "sex", attr_actors, scaling = "std") +
        remulate::difference(0.4, "age", attr_actors, scaling = "std") +
        remulate::outdegreeReceiver(0.5, scaling = "std") +
        remulate::inertia(0.4, scaling = "std")
        
    } else if (scenario_name == "HARD") {
      effects1 <- ~
        remulate::baseline(-6) +
        remulate::difference(0.3, "sex", attr_actors, scaling = "std") +
        remulate::difference(0.4, "age", attr_actors, scaling = "std") +
        remulate::outdegreeReceiver(0.3, scaling = "std") +
        remulate::inertia(0.25, scaling = "std")
      
      effects2 <- ~
        remulate::baseline(-5) +
        remulate::difference(0.3, "sex", attr_actors, scaling = "std") +
        remulate::difference(0.4, "age", attr_actors, scaling = "std") +
        remulate::outdegreeReceiver(0.4, scaling = "std") +
        remulate::inertia(0.35, scaling = "std")
        
    } else if (scenario_name == "EHARD") {
      effects1 <- ~
        remulate::baseline(-6) +
        remulate::difference(0.3, "sex", attr_actors, scaling = "std") +
        remulate::difference(0.4, "age", attr_actors, scaling = "std") +
        remulate::outdegreeReceiver(0.2, scaling = "std") +
        remulate::inertia(0.25, scaling = "std")
      
      effects2 <- ~
        remulate::baseline(-5) +
        remulate::difference(0.3, "sex", attr_actors, scaling = "std") +
        remulate::difference(0.4, "age", attr_actors, scaling = "std") +
        remulate::outdegreeReceiver(0.45, scaling = "std") +
        remulate::inertia(0.4, scaling = "std")
    }
    
    # Initialize and simulate events
    initialREH <- data.frame(time = 1, sender = 1, receiver = 1)
    all_events <- NULL
    actual_states <- c()  # Track the actual states for each event
    
    for (i in seq_along(state_lengths)) {
      n_segment_events <- state_lengths[i]
      current_state <- state_values[i]
      current_effects <- if (current_state == 1) effects1 else effects2
      
      # Skip very short segments (less than 2 events) to avoid edge cases
      if (n_segment_events < 2) {
        next
      }
      
      sim <- tryCatch({
        remulate::remulateTie(
          effects = current_effects,
          actors = actors,
          events = n_segment_events,
          endTime = 100000,
          initial = initialREH
        )
      }, error = function(e) {
        message("    Warning: Skipping segment due to error: ", e$message)
        NULL
      })
      
      if (!is.null(sim) && nrow(sim) > 0) {
        all_events <- rbind(all_events, sim)
        actual_states <- c(actual_states, rep(current_state, nrow(sim)))
        initialREH <- all_events
      }
    }
    
    # Check if we have enough events
    if (is.null(all_events) || nrow(all_events) < 100) {
      message("    Warning: Not enough events generated, retrying iteration...")
      next
    }
    
    # Prepare event data frame
    events_df <- as.data.frame(all_events)
    events_df$sender <- as.character(events_df$sender)
    events_df$receiver <- as.character(events_df$receiver)
    events_df$true_state <- actual_states[seq_len(nrow(events_df))]
    
    # Compute time differences
    time_differences <- diff(events_df[[1]])
    adjusted_differences <- time_differences
    last_nonzero_diff <- NA
    
    for (j in 1:length(time_differences)) {
      if (time_differences[j] == 0) {
        adjusted_differences[j] <- last_nonzero_diff
      } else {
        last_nonzero_diff <- time_differences[j]
      }
    }
    
    adjusted_differences[is.na(adjusted_differences)] <- 0
    adjusted_differences0 <- c(0, adjusted_differences)
    events_df$Timedifferencees <- adjusted_differences0
    events_df$Timedifferencees[1] <- events_df$time[1]
    
    # Prepare HMM data
    dt <- events_df$Timedifferencees
    HMMdf <- data.frame(ID = rep("a", length(dt)), step = dt, angle = NA)
    hmm_data <- prepData(data = HMMdf, coordNames = NULL)
    
    # Storage for results
    rep_results <- list(iter = iter, n_events = nrow(events_df), 
                        hmm_results = list(), rem_results = list())
    
    # Fit HMM with different numbers of states
    for (m in states_to_fit) {
      cat("    Fitting HMM with", m, "state(s)... ")
      
      hmm_fit <- fit_hmm_safe(hmm_data, m)
      
      if (!is.null(hmm_fit)) {
        bic <- AIC(hmm_fit, k = log(nrow(hmm_data)))
        aic <- AIC(hmm_fit)
        
        # Handle 1-state case differently
        if (m == 1) {
          conf_metrics <- list(avg_confidence = 1.0, pct_uncertain_tau1 = 0, pct_uncertain_tau2 = 0)
          viterbi_states <- rep(1, nrow(events_df))
          accuracy <- NA
        } else {
          conf_metrics <- compute_confidence_metrics(hmm_fit, tau1, tau2)
          viterbi_states <- momentuHMM::viterbi(hmm_fit)
          accuracy <- calculate_accuracy(events_df$true_state, viterbi_states, m)
        }
        
        rep_results$hmm_results[[paste0("m", m)]] <- list(
          n_states = m, bic = bic, aic = aic,
          avg_confidence = conf_metrics$avg_confidence,
          pct_uncertain_tau1 = conf_metrics$pct_uncertain_tau1,
          pct_uncertain_tau2 = conf_metrics$pct_uncertain_tau2,
          accuracy = accuracy,
          step_params = hmm_fit$mle$step,
          gamma = if (m > 1) hmm_fit$mle$gamma else NULL
        )
        
        cat("Done (BIC:", round(bic, 2), ")\n")
        
        # For 2-state model, fit REM models
        if (m == 2) {
          events_df$predicted_state <- viterbi_states
          means_by_group <- aggregate(Timedifferencees ~ predicted_state, 
                                      data = events_df, FUN = mean)
          sorted_index <- order(means_by_group$Timedifferencees, decreasing = TRUE)
          events_df$state_indicator <- ifelse(events_df$predicted_state == sorted_index[2], 1, 0)
          
          reh_tie <- remify::remify(edgelist = events_df, model = "tie", 
                                    actors = attr_actors$name, directed = TRUE, origin = 0)
          
          # Model 1: Simple REM
          stats1 <- ~ 1 + difference("sex", scaling = "std") +
            difference("age", scaling = "std") +
            outdegreeReceiver(scaling = "std") +
            inertia(scaling = "std")
          
          out1 <- remstats(reh = reh_tie, tie_effects = stats1, attr_actors = attr_actors)
          fit1 <- remstimate::remstimate(reh = reh_tie, stats = out1, method = "MLE")
          summary1 <- summary(fit1)
          
          rep_results$rem_results[["simple"]] <- list(
            model = "Simple REM", coefs = summary1$coefsTab,
            bic = summary1$BIC, aic = summary1$AIC
          )
          
          # Model 2: REM with State effect
          stats2 <- ~ 1 + difference("sex", scaling = "std") +
            difference("age", scaling = "std") +
            outdegreeReceiver(scaling = "std") +
            inertia(scaling = "std") +
            (event(x = events_df$state_indicator, "PredictedState1"))
          
          out2 <- remstats(reh = reh_tie, tie_effects = stats2, attr_actors = attr_actors)
          fit2 <- remstimate::remstimate(reh = reh_tie, stats = out2, method = "MLE")
          summary2 <- summary(fit2)
          
          rep_results$rem_results[["state"]] <- list(
            model = "REM + State", coefs = summary2$coefsTab,
            bic = summary2$BIC, aic = summary2$AIC
          )
          
          # Model 3: REM with Interaction effect
          stats3 <- ~ 1 + difference("sex", scaling = "std") +
            difference("age", scaling = "std") +
            (outdegreeReceiver(scaling = "std") + inertia(scaling = "std")) :
            (event(x = events_df$state_indicator, "PredictedState1"))
          
          out3 <- remstats(reh = reh_tie, tie_effects = stats3, attr_actors = attr_actors)
          fit3 <- remstimate::remstimate(reh = reh_tie, stats = out3, method = "MLE")
          summary3 <- summary(fit3)
          
          rep_results$rem_results[["interaction"]] <- list(
            model = "REM + Interaction", coefs = summary3$coefsTab,
            bic = summary3$BIC, aic = summary3$AIC
          )
        }
      } else {
        cat("Failed\n")
        rep_results$hmm_results[[paste0("m", m)]] <- list(
          n_states = m, bic = NA, aic = NA, avg_confidence = NA, 
          pct_uncertain = NA, accuracy = NA
        )
      }
    }
    
    replications[[iter]] <- rep_results
  }
  
  all_results[[scenario_name]] <- list(scenario = scenario_name, replications = replications)
  all_hmm_agg[[scenario_name]] <- aggregate_hmm_results(scenario_name, replications)
  all_rem_agg[[scenario_name]] <- aggregate_rem_results(replications)
}

# -----------------------------------------------------------------------------
# 7. Results Summary
# -----------------------------------------------------------------------------

cat("\n\n", strrep("=", 70), "\n")
cat("  AGGREGATED RESULTS\n")
cat(strrep("=", 70), "\n")

combined_hmm <- do.call(rbind, all_hmm_agg)

cat("\n--- HMM Model Selection Summary ---\n\n")
print(knitr::kable(combined_hmm[, c("Scenario", "States", "BIC_mean", "BIC_sd", 
                                     "Confidence_mean", "Uncertainty_tau1_mean", 
                                     "Uncertainty_tau2_mean", "Accuracy_mean")],
                   digits = 2,
                   col.names = c("Scenario", "States", "BIC (mean)", "BIC (sd)",
                                 "Confidence", paste0("% Unc (tau=", tau1, ")"), 
                                 paste0("% Unc (tau=", tau2, ")"), "Accuracy")))

cat("\n\n--- REM Coefficient Estimates (2-State Model) ---\n")

for (scenario_name in names(all_rem_agg)) {
  cat("\n", strrep("-", 50), "\n")
  cat("Scenario:", scenario_name, "\n")
  cat(strrep("-", 50), "\n")
  
  rem_results <- all_rem_agg[[scenario_name]]
  
  if (!is.null(rem_results)) {
    for (model_name in names(rem_results)) {
      cat("\n  Model:", rem_results[[model_name]]$model, "\n")
      cat("  BIC: ", round(rem_results[[model_name]]$bic_mean, 2), 
          " (", round(rem_results[[model_name]]$bic_sd, 2), ")\n", sep = "")
      cat("\n")
      print(knitr::kable(rem_results[[model_name]]$coefficients, digits = 3))
    }
  }
}

# -----------------------------------------------------------------------------
# 8. Generate Figures
# -----------------------------------------------------------------------------

cat("\n\n--- Generating Figures ---\n")

p_bic <- create_bic_comparison_plot(all_hmm_agg, R)
print(p_bic)
ggsave("HMMREM_scenario_comparison.png", p_bic, width = 10, height = 6, dpi = 300)
cat("Saved: HMMREM_scenario_comparison.png\n")

p_conf <- create_confidence_plot(all_hmm_agg, tau1, tau2)
print(p_conf)
ggsave("HMMREM_confidence_comparison.png", p_conf, width = 8, height = 6, dpi = 300)
cat("Saved: HMMREM_confidence_comparison.png\n")

# -----------------------------------------------------------------------------
# 9. Final Summary Table
# -----------------------------------------------------------------------------

cat("\n\n", strrep("=", 70), "\n")
cat("  FINAL SUMMARY TABLE\n")
cat(strrep("=", 70), "\n\n")

summary_table <- combined_hmm %>%
  filter(States == 2) %>%
  select(Scenario, BIC_mean, BIC_sd, Confidence_mean, Confidence_sd, 
         Uncertainty_tau1_mean, Uncertainty_tau2_mean, Accuracy_mean, Accuracy_sd) %>%
  mutate(
    BIC = paste0(round(BIC_mean, 1), " (", round(BIC_sd, 1), ")"),
    Confidence = paste0(round(Confidence_mean * 100, 1), "% (", 
                        round(Confidence_sd * 100, 1), ")"),
    Unc_tau1 = paste0(round(Uncertainty_tau1_mean, 1), "%"),
    Unc_tau2 = paste0(round(Uncertainty_tau2_mean, 1), "%"),
    Accuracy = paste0(round(Accuracy_mean, 1), "% (", round(Accuracy_sd, 1), ")")
  ) %>%
  select(Scenario, BIC, Confidence, Unc_tau1, Unc_tau2, Accuracy)

cat("2-State HMM Performance Summary:\n\n")
print(knitr::kable(summary_table, align = "lcccc",
                   col.names = c("Scenario", "BIC", "Confidence", 
                                 paste0("Unc (tau=", tau1, ")"),
                                 paste0("Unc (tau=", tau2, ")"), "Accuracy")))

cat("\n\nOptimal Number of States (by BIC):\n")
for (scenario_name in names(all_hmm_agg)) {
  best_row <- all_hmm_agg[[scenario_name]][which.min(all_hmm_agg[[scenario_name]]$BIC_mean), ]
  cat("  ", scenario_name, ": ", best_row$States, " state(s) (BIC = ", 
      round(best_row$BIC_mean, 1), ")\n", sep = "")
}

cat("\n\nAnalysis complete!\n")

# -----------------------------------------------------------------------------
# 10. Save results to global environment for Excel export
# -----------------------------------------------------------------------------

# Ensure all results are available in global environment
assign("all_hmm_agg", all_hmm_agg, envir = .GlobalEnv)
assign("all_rem_agg", all_rem_agg, envir = .GlobalEnv)
assign("all_results", all_results, envir = .GlobalEnv)
assign("scenarios", scenarios, envir = .GlobalEnv)
assign("R", R, envir = .GlobalEnv)
assign("n_actors", n_actors, envir = .GlobalEnv)
assign("states_to_fit", states_to_fit, envir = .GlobalEnv)
assign("tau1", tau1, envir = .GlobalEnv)
assign("tau2", tau2, envir = .GlobalEnv)
