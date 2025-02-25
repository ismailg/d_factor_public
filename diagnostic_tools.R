###############################################################################
# Parameter Identifiability Diagnostic Tools
# These functions help analyze which parameters can be reliably recovered
###############################################################################

library(ggplot2)
library(gridExtra)
library(reshape2)
library(viridis)

###############################################################################
# 1) Parameter Slices Analysis
###############################################################################

# Generate model predictions along parameter slices
# What it does: Varies one parameter at a time around its estimated value while holding others fixed, calculating the negative log-likelihood (NLL).
# Purpose: Identifies parameters with "flat" likelihood surfaces (hard to estimate) versus those with clear minima (well-identified).
# Visualization: Line plots showing NLL changes with dashed (original) and dotted (optimal) vertical lines.
generate_parameter_slices <- function(best_params, data, step_size=10) {
  all_params <- c("alpha", "temp", "envy", "guilt")
  slice_results <- list()
  
  for (param in all_params) {
    current_value <- best_params[[param]]
    
    # Skip if parameter is zero (likely fixed)
    if (is.null(current_value) || current_value == 0) next
    
    # Create range around the parameter value
    if (param == "alpha") {
      param_range <- seq(max(0.01, current_value/3), min(0.99, current_value*3), length.out=step_size)
    } else if (param == "temp") {
      param_range <- seq(max(0.1, current_value/3), current_value*3, length.out=step_size)
    } else if (param == "envy") {
      param_range <- seq(max(0, current_value/3), current_value*3, length.out=step_size)
    } else if (param == "guilt") {
      param_range <- seq(max(0, current_value/3), current_value*3, length.out=step_size)
    }
    
    # Calculate NLL across the range
    nll_values <- numeric(length(param_range))
    
    for (i in seq_along(param_range)) {
      # Make a copy of the best parameters
      modified_params <- best_params
      modified_params[[param]] <- param_range[i]
      
      # Extract free parameters for model evaluation
      params_free <- unlist(modified_params[all_params])
      
      # Calculate NLL
      nll_values[i] <- direct_learning_model(params_free, data, list())
    }
    
    # Store results
    slice_results[[param]] <- list(
      parameter = param,
      values = param_range,
      nll = nll_values,
      best_value = current_value,
      min_value = param_range[which.min(nll_values)]
    )
  }
  
  return(slice_results)
}

# Plot parameter slices
plot_parameter_slices <- function(slice_results) {
  if (length(slice_results) == 0) {
    message("No slice results to plot")
    return(NULL)
  }
  
  plots <- list()
  
  for (param in names(slice_results)) {
    result <- slice_results[[param]]
    
    # Normalize NLL values to start at 0 for better visualization
    nll_norm <- result$nll - min(result$nll)
    
    # Create data frame for plotting
    df <- data.frame(
      value = result$values,
      nll = nll_norm
    )
    
    # Create plot
    p <- ggplot(df, aes(x = value, y = nll)) +
      geom_line(size = 1) +
      geom_point() +
      geom_vline(xintercept = result$best_value, linetype = "dashed", color = "red") +
      geom_vline(xintercept = result$min_value, linetype = "dotted", color = "blue") +
      labs(title = paste("NLL Profile for", param),
           x = param,
           y = "Delta NLL") +
      theme_minimal()
    
    plots[[param]] <- p
  }
  
  # Arrange plots in a grid
  if (length(plots) > 0) {
    grid_size <- length(plots)
    n_rows <- ceiling(sqrt(grid_size))
    n_cols <- ceiling(grid_size / n_rows)
    
    grid.arrange(grobs = plots, ncol = n_cols)
  }
}

###############################################################################
# 2) Parameter Interaction Analysis
###############################################################################
# What it does: Creates 2D grids of parameter pairs to analyze how NLL changes when two parameters vary simultaneously.
# Purpose: Reveals correlated parameters that might be non-identifiable when considered together.
# Visualization: Heatmaps with a red dot marking the original parameter values.

# Analyze pairwise parameter interactions
analyze_parameter_interactions <- function(best_params, data, n_steps=10) {
  all_params <- c("alpha", "temp", "envy", "guilt")
  
  # Only analyze parameters that are not fixed
  free_params <- all_params[!sapply(best_params[all_params], is.null)]
  
  if (length(free_params) < 2) {
    message("Fewer than 2 free parameters. No interactions to analyze.")
    return(NULL)
  }
  
  # Generate all pairwise combinations
  param_pairs <- combn(free_params, 2, simplify = FALSE)
  interaction_results <- list()
  
  # Create parameter grid for each pair
  for (pair in param_pairs) {
    param1 <- pair[1]
    param2 <- pair[2]
    
    # Center values around best estimates
    val1 <- best_params[[param1]]
    val2 <- best_params[[param2]]
    
    # Define ranges
    range1 <- get_parameter_range(param1, val1, n_steps)
    range2 <- get_parameter_range(param2, val2, n_steps)
    
    # Create grid
    nll_grid <- matrix(NA, nrow = length(range1), ncol = length(range2))
    
    # Calculate NLL for each grid point
    for (i in seq_along(range1)) {
      for (j in seq_along(range2)) {
        # Make a copy of the best parameters
        modified_params <- best_params
        modified_params[[param1]] <- range1[i]
        modified_params[[param2]] <- range2[j]
        
        # Calculate NLL
        nll_grid[i, j] <- direct_learning_model(
          unlist(modified_params[all_params]), data, list())
      }
    }
    
    # Normalize NLL values to start at 0
    nll_grid_norm <- nll_grid - min(nll_grid, na.rm = TRUE)
    
    # Store results
    interaction_results[[paste(param1, param2, sep = "_")]] <- list(
      param1 = param1,
      param2 = param2,
      range1 = range1,
      range2 = range2,
      nll_grid = nll_grid_norm,
      best_value1a = val1,
      best_value2 = val2
    )
  }
  
  return(interaction_results)
}

# Helper function to define parameter ranges
get_parameter_range <- function(param, value, n_steps) {
  if (param == "alpha") {
    return(seq(max(0.01, value/2), min(0.99, value*2), length.out = n_steps))
  } else if (param == "temp") {
    return(seq(max(0.1, value/2), value*2, length.out = n_steps))
  } else if (param == "envy") {
    return(seq(max(0, value/2), value*2, length.out = n_steps))
  } else if (param == "guilt") {
    return(seq(max(0, value/2), value*2, length.out = n_steps))
  }
}

# Plot parameter interactions
plot_parameter_interactions <- function(interaction_results) {
  if (length(interaction_results) == 0) {
    message("No interaction results to plot")
    return(NULL)
  }
  
  plots <- list()
  
  for (interaction_name in names(interaction_results)) {
    result <- interaction_results[[interaction_name]]
    
    # Convert matrix to data frame for plotting
    df <- expand.grid(
      x = result$range1,
      y = result$range2
    )
    df$z <- as.vector(result$nll_grid)
    
    # Create plot
    p <- ggplot(df, aes(x = x, y = y, fill = z)) +
      geom_tile() +
      scale_fill_viridis_c(name = "Delta NLL", 
                           limits = c(0, min(max(df$z, na.rm = TRUE), 10))) +
      geom_point(aes(x = result$best_value1, y = result$best_value2), 
                 color = "red", size = 3) +
      labs(title = paste("Parameter Interaction:", 
                         result$param1, "vs", result$param2),
           x = result$param1,
           y = result$param2) +
      theme_minimal()
    
    plots[[interaction_name]] <- p
  }
  
  # Arrange plots in a grid
  if (length(plots) > 0) {
    grid_size <- length(plots)
    n_rows <- ceiling(sqrt(grid_size))
    n_cols <- ceiling(grid_size / n_rows)
    
    grid.arrange(grobs = plots, ncol = n_cols)
  }
}

###############################################################################
# 3) Parameter Recovery Matrix
###############################################################################
# What it does: Simulates data with known parameters and checks how well the model can recover them.
# Metrics: Correlation between true and recovered parameters
# Mean Absolute Error (MAE) in recovery
# Visualization: Heatmap of correlations and bar plot of MAEs.

# Generate a parameter recovery quality matrix
generate_recovery_matrix <- function(recovery_results) {
  if (is.null(recovery_results$summary_df) || nrow(recovery_results$summary_df) < 3) {
    message("Not enough recovery results to generate matrix")
    return(NULL)
  }
  
  all_params <- c("alpha", "temp", "envy", "guilt")
  df <- recovery_results$summary_df
  
  # Calculate correlations
  correlations <- matrix(NA, nrow = length(all_params), ncol = length(all_params))
  rownames(correlations) <- all_params
  colnames(correlations) <- all_params
  
  for (i in seq_along(all_params)) {
    for (j in seq_along(all_params)) {
      true_col <- paste0("true_", all_params[i])
      fit_col <- paste0("fit_", all_params[j])
      
      if (true_col %in% colnames(df) && fit_col %in% colnames(df)) {
        valid_idx <- !is.na(df[[true_col]]) & !is.na(df[[fit_col]])
        
        if (sum(valid_idx) >= 3) {
          correlations[i, j] <- cor(df[[true_col]][valid_idx], 
                                    df[[fit_col]][valid_idx])
        }
      }
    }
  }
  
  # Calculate mean absolute errors
  mae <- matrix(NA, nrow = length(all_params), ncol = 1)
  rownames(mae) <- all_params
  colnames(mae) <- "MAE"
  
  for (i in seq_along(all_params)) {
    param <- all_params[i]
    error_col <- paste0(param, "_error")
    
    if (error_col %in% colnames(df)) {
      mae[i, 1] <- mean(df[[error_col]], na.rm = TRUE)
    }
  }
  
  # Return both matrices
  list(
    correlations = correlations,
    mae = mae
  )
}

# Plot parameter recovery matrix
plot_recovery_matrix <- function(recovery_matrix) {
  if (is.null(recovery_matrix) || is.null(recovery_matrix$correlations)) {
    message("No recovery matrix to plot")
    return(NULL)
  }
  
  # Convert correlation matrix to long format
  corr_matrix <- recovery_matrix$correlations
  melted_corr <- melt(corr_matrix)
  names(melted_corr) <- c("True", "Estimated", "Correlation")
  
  # Create correlation heatmap
  p1 <- ggplot(melted_corr, aes(x = True, y = Estimated, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-1, 1)) +
    geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
    labs(title = "Parameter Recovery Correlation Matrix",
         x = "True Parameter",
         y = "Estimated Parameter") +
    theme_minimal()
  
  # Convert MAE to long format
  mae_matrix <- recovery_matrix$mae
  melted_mae <- melt(mae_matrix)
  names(melted_mae) <- c("Parameter", "Metric", "Value")
  
  # Create MAE bar plot
  p2 <- ggplot(melted_mae, aes(x = Parameter, y = Value)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    geom_text(aes(label = round(Value, 3)), vjust = -0.5) +
    labs(title = "Mean Absolute Error by Parameter",
         x = "Parameter",
         y = "MAE") +
    theme_minimal()
  
  # Arrange plots
  grid.arrange(p1, p2, ncol = 1, heights = c(2, 1))
}

###############################################################################
# 4) Posterior Predictive Checks
###############################################################################

# Generate posterior predictive checks
generate_posterior_checks <- function(data, best_params, n_simulations=10) {
  # Extract parameters
  all_params <- c("alpha", "temp", "envy", "guilt")
  
  # Run simulations with best parameters
  simulations <- list()
  
  for (i in 1:n_simulations) {
    # Set seed for reproducibility but different for each simulation
    set.seed(i)
    
    # Run simulation
    sim_data <- simulate_2games_data_mf_hmm(
      best_params, 
      game_size = max(data$roundNum, na.rm = TRUE), 
      playerId = i
    )
    
    simulations[[i]] <- sim_data
  }
  
  # Extract key statistics from real data
  real_stats <- extract_statistics(data)
  
  # Extract statistics from simulations
  sim_stats <- lapply(simulations, extract_statistics)
  
  # Combine all simulation statistics
  sim_stats_combined <- list()
  for (stat_name in names(real_stats)) {
    sim_values <- sapply(sim_stats, function(s) s[[stat_name]])
    sim_stats_combined[[stat_name]] <- sim_values
  }
  
  # Return all statistics
  list(
    real_stats = real_stats,
    sim_stats = sim_stats_combined,
    simulations = simulations
  )
}

# Helper function to extract key statistics from data
extract_statistics <- function(data) {
  stats <- list()
  
  # Calculate average return proportion by game
  if ("gameNum.f" %in% names(data) && "investment" %in% names(data) && "return" %in% names(data)) {
    game_stats <- aggregate(
      cbind(return, investment) ~ gameNum.f, 
      data = data, 
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    
    game_stats$return_prop <- game_stats$return / (3 * game_stats$investment)
    stats$mean_return_prop_by_game <- game_stats$return_prop
  }
  
  # Calculate return proportion by investor state if available
  if ("investor_state" %in% names(data) && "investment" %in% names(data) && "return" %in% names(data)) {
    state_stats <- aggregate(
      cbind(return, investment) ~ investor_state, 
      data = data, 
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    
    state_stats$return_prop <- state_stats$return / (3 * state_stats$investment)
    stats$mean_return_prop_by_state <- state_stats$return_prop
    names(stats$mean_return_prop_by_state) <- state_stats$investor_state
  }
  
  # Calculate average return proportion by round
  if ("roundNum" %in% names(data) && "investment" %in% names(data) && "return" %in% names(data)) {
    round_stats <- aggregate(
      cbind(return, investment) ~ roundNum, 
      data = data, 
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    
    round_stats$return_prop <- round_stats$return / (3 * round_stats$investment)
    stats$mean_return_prop_by_round <- round_stats$return_prop
  }
  
  return(stats)
}

# Plot posterior predictive checks
plot_posterior_checks <- function(posterior_checks) {
  if (is.null(posterior_checks) || is.null(posterior_checks$real_stats)) {
    message("No posterior check results to plot")
    return(NULL)
  }
  
  plots <- list()
  
  # Plot return proportion by game
  if (!is.null(posterior_checks$real_stats$mean_return_prop_by_game)) {
    real_values <- posterior_checks$real_stats$mean_return_prop_by_game
    sim_values <- posterior_checks$sim_stats$mean_return_prop_by_game
    
    # Create data frame for plotting
    n_games <- length(real_values)
    n_sims <- ncol(sim_values)
    
    df_sim <- data.frame(
      game = rep(1:n_games, each = n_sims),
      return_prop = as.vector(sim_values),
      type = "Simulated"
    )
    
    df_real <- data.frame(
      game = 1:n_games,
      return_prop = real_values,
      type = "Real"
    )
    
    df <- rbind(df_sim, df_real)
    
    # Create plot
    p1 <- ggplot(df, aes(x = factor(game), y = return_prop, color = type)) +
      geom_boxplot(data = subset(df, type == "Simulated")) +
      geom_point(data = subset(df, type == "Real"), size = 3) +
      labs(title = "Return Proportion by Game",
           x = "Game",
           y = "Return Proportion") +
      theme_minimal()
    
    plots[["by_game"]] <- p1
  }
  
  # Plot return proportion by round
  if (!is.null(posterior_checks$real_stats$mean_return_prop_by_round)) {
    real_values <- posterior_checks$real_stats$mean_return_prop_by_round
    sim_values <- posterior_checks$sim_stats$mean_return_prop_by_round
    
    # Create data frame for plotting
    n_rounds <- length(real_values)
    n_sims <- ncol(sim_values)
    
    df_sim <- data.frame(
      round = rep(1:n_rounds, each = n_sims),
      return_prop = as.vector(sim_values),
      type = "Simulated"
    )
    
    df_real <- data.frame(
      round = 1:n_rounds,
      return_prop = real_values,
      type = "Real"
    )
    
    df <- rbind(df_sim, df_real)
    
    # Create plot
    p2 <- ggplot(df, aes(x = factor(round), y = return_prop, color = type)) +
      geom_boxplot(data = subset(df, type == "Simulated")) +
      geom_point(data = subset(df, type == "Real"), size = 3) +
      labs(title = "Return Proportion by Round",
           x = "Round",
           y = "Return Proportion") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    plots[["by_round"]] <- p2
  }
  
  # Plot return proportion by investor state
  if (!is.null(posterior_checks$real_stats$mean_return_prop_by_state)) {
    real_values <- posterior_checks$real_stats$mean_return_prop_by_state
    sim_values <- posterior_checks$sim_stats$mean_return_prop_by_state
    
    if (!is.null(sim_values)) {
      # Create data frame for plotting
      n_states <- length(real_values)
      n_sims <- ncol(sim_values)
      
      df_sim <- data.frame(
        state = rep(names(real_values), each = n_sims),
        return_prop = as.vector(sim_values),
        type = "Simulated"
      )
      
      df_real <- data.frame(
        state = names(real_values),
        return_prop = real_values,
        type = "Real"
      )
      
      df <- rbind(df_sim, df_real)
      
      # Create plot
      p3 <- ggplot(df, aes(x = state, y = return_prop, color = type)) +
        geom_boxplot(data = subset(df, type == "Simulated")) +
        geom_point(data = subset(df, type == "Real"), size = 3) +
        labs(title = "Return Proportion by Investor State",
             x = "Investor State",
             y = "Return Proportion") +
        theme_minimal()
      
      plots[["by_state"]] <- p3
    }
  }
  
  # Arrange plots in a grid
  if (length(plots) > 0) {
    grid_size <- length(plots)
    n_rows <- grid_size
    
    grid.arrange(grobs = plots, ncol = 1)
  }
}

###############################################################################
# 5) Comprehensive Diagnostics Suite
###############################################################################

# Run a comprehensive set of diagnostics on a dataset and model fit
run_comprehensive_diagnostics <- function(data, fit_result, n_recovery=10, verbose=FALSE) {
  if (verbose) cat("Starting comprehensive diagnostics...\n")
  
  # Extract best parameters
  best_params <- fit_result$params
  
  # 1. Parameter Slices Analysis
  if (verbose) cat("Generating parameter slices...\n")
  slice_results <- generate_parameter_slices(best_params, data)
  
  # 2. Parameter Interaction Analysis
  if (verbose) cat("Analyzing parameter interactions...\n")
  interaction_results <- analyze_parameter_interactions(best_params, data)
  
  # 3. Parameter Recovery Analysis
  if (verbose) cat("Running parameter recovery tests...\n")
  recovery_results <- param_recovery_experiment(
    n_sims = n_recovery, 
    game_size = max(data$roundNum, na.rm = TRUE),
    param_fixed = list(),
    verbose = verbose
  )
  recovery_matrix <- generate_recovery_matrix(recovery_results)
  
  # 4. Posterior Predictive Checks
  if (verbose) cat("Generating posterior predictive checks...\n")
  posterior_checks <- generate_posterior_checks(data, best_params)
  
  # 5. Generate Diagnostic Plots
  if (verbose) cat("Generating diagnostic plots...\n")
  
  # Parameter slices
  if (verbose) cat("  - Parameter slices plot\n")
  plot_parameter_slices(slice_results)
  
  # Parameter interactions
  if (verbose) cat("  - Parameter interactions plot\n")
  plot_parameter_interactions(interaction_results)
  
  # Recovery matrix
  if (verbose) cat("  - Parameter recovery matrix\n")
  plot_recovery_matrix(recovery_matrix)
  
  # Posterior checks
  if (verbose) cat("  - Posterior predictive checks\n")
  plot_posterior_checks(posterior_checks)
  
  # 6. Generate Recommendations
  if (verbose) cat("Generating parameter recommendations...\n")
  
  # Analyze which parameters are most identifiable
  param_quality <- list()
  
  # From slice results: assess curvature
  if (!is.null(slice_results)) {
    for (param in names(slice_results)) {
      result <- slice_results[[param]]
      # Rough measure of curvature - variance in NLL values
      curvature <- var(result$nll)
      param_quality[[param]] <- list(
        curvature = curvature
      )
    }
  }
  
  # From recovery results: get correlations
  if (!is.null(recovery_matrix)) {
    for (param in rownames(recovery_matrix$correlations)) {
      # Diagonal element is correlation between true and estimated
      corr <- recovery_matrix$correlations[param, param]
      mae <- recovery_matrix$mae[param, 1]
      
      if (!is.null(param_quality[[param]])) {
        param_quality[[param]]$correlation <- corr
        param_quality[[param]]$mae <- mae
      } else {
        param_quality[[param]] <- list(
          correlation = corr,
          mae = mae,
          curvature = NA
        )
      }
    }
  }
  
  # Calculate overall quality score
  quality_scores <- numeric(length(param_quality))
  names(quality_scores) <- names(param_quality)
  
  for (param in names(param_quality)) {
    # Higher correlation, higher curvature, lower MAE is better
    correlation <- param_quality[[param]]$correlation
    curvature <- param_quality[[param]]$curvature
    mae <- param_quality[[param]]$mae
    
    # Normalize all to [0,1] range
    corr_score <- if (!is.na(correlation)) max(0, min(1, (correlation + 1) / 2)) else 0
    curv_score <- if (!is.na(curvature) && !is.null(curvature)) {
      all_curvatures <- sapply(param_quality, function(x) x$curvature)
      all_curvatures <- all_curvatures[!is.na(all_curvatures) & !is.null(all_curvatures)]
      if (length(all_curvatures) > 0) {
        min_curv <- min(all_curvatures)
        max_curv <- max(all_curvatures)
        if (max_curv > min_curv) {
          (curvature - min_curv) / (max_curv - min_curv)
        } else {
          0.5  # Default if all curvatures are the same
        }
      } else {
        0.5  # Default if no valid curvatures
      }
    } else {
      0.5  # Default if curvature not available
    }
    
    mae_score <- if (!is.na(mae) && !is.null(mae)) {
      all_maes <- sapply(param_quality, function(x) x$mae)
      all_maes <- all_maes[!is.na(all_maes) & !is.null(all_maes)]
      if (length(all_maes) > 0) {
        min_mae <- min(all_maes)
        max_mae <- max(all_maes)
        if (max_mae > min_mae) {
          1 - (mae - min_mae) / (max_mae - min_mae)  # Invert so lower MAE = higher score
        } else {
          0.5  # Default if all MAEs are the same
        }
      } else {
        0.5  # Default if no valid MAEs
      }
    } else {
      0.5  # Default if MAE not available
    }
    
    # Weighted average of scores
    quality_scores[param] <- 0.5*corr_score + 0.25*curv_score + 0.25*mae_score
  }
  
  # Generate recommendations
  recommendations <- list()
  
  for (param in names(quality_scores)) {
    score <- quality_scores[param]
    
    if (score >= 0.7) {
      recommendation <- "KEEP: This parameter is well-identified and should be kept free."
    } else if (score >= 0.4) {
      recommendation <- "CONSIDER: This parameter is moderately identified. Consider group-level constraints or informative priors."
    } else {
      recommendation <- "FIX: This parameter is poorly identified. Fix to a sensible value from literature."
      
      # Suggest value
      if (param == "alpha") {
        value <- 0.3  # Common learning rate value
      } else if (param == "temp") {
        value <- 1.5  # Common temperature value
      } else if (param == "envy") {
        value <- 1.0  # Common envy parameter
      } else if (param == "guilt") {
        value <- 0.5  # Common guilt parameter
      } else {
        value <- best_params[[param]]  # Default to current value
      }
      
      recommendation <- paste0(recommendation, " Suggested value: ", value)
    }
    
    recommendations[[param]] <- list(
      score = score,
      recommendation = recommendation,
      current_value = best_params[[param]]
    )
  }
  
  # 7. Print summary report
  if (verbose) {
    cat("\nDiagnostic Summary:\n")
    cat("------------------\n")
    
    for (param in names(recommendations)) {
      cat(param, ":\n")
      cat("  Identifiability Score:", round(recommendations[[param]]$score, 2), "\n")
      cat("  Current Value:", round(recommendations[[param]]$current_value, 3), "\n")
      cat("  Recommendation:", recommendations[[param]]$recommendation, "\n\n")
    }
  }
  
  # Return all results
  list(
    slice_results = slice_results,
    interaction_results = interaction_results,
    recovery_results = recovery_results,
    recovery_matrix = recovery_matrix,
    posterior_checks = posterior_checks,
    quality_scores = quality_scores,
    recommendations = recommendations
  )
}

###############################################################################
# Main function to combine all analyses
###############################################################################

analyze_model_identifiability <- function(data, n_recovery=10, verbose=TRUE) {
  # Step 1: Fit the model with all parameters free
  if (verbose) cat("Fitting model with all parameters free...\n")
  fit_all_free <- fit_direct_model_robust(data, param_fixed=list(), 
                                          n_multistart=10, verbose=verbose)
  
  # Step 2: Run comprehensive diagnostics
  if (verbose) cat("Running comprehensive diagnostics...\n")
  diagnostics <- run_comprehensive_diagnostics(data, fit_all_free, 
                                               n_recovery=n_recovery, 
                                               verbose=verbose)
  
  # Step 3: Get recommendations for which parameters to fix
  recommendations <- diagnostics$recommendations
  
  # Find parameters to fix (score < 0.4)
  params_to_fix <- list()
  for (param in names(recommendations)) {
    if (recommendations[[param]]$score < 0.4) {
      # Extract suggested value
      rec <- recommendations[[param]]$recommendation
      value_str <- sub(".*Suggested value: ([0-9.]+).*", "\\1", rec)
      value <- as.numeric(value_str)
      
      params_to_fix[[param]] <- value
    }
  }
  
  # Step 4: Refit with fixed parameters
  if (length(params_to_fix) > 0) {
    if (verbose) {
      cat("\nRefitting model with fixed parameters:\n")
      for (param in names(params_to_fix)) {
        cat("  ", param, "=", params_to_fix[[param]], "\n")
      }
    }
    
    fit_final <- fit_direct_model_robust(data, param_fixed=params_to_fix, 
                                         n_multistart=10, verbose=verbose)
  } else {
    if (verbose) cat("\nNo parameters need to be fixed. Using original fit.\n")
    fit_final <- fit_all_free
  }
  
  # Step 5: Compare models if parameters were fixed
  if (length(params_to_fix) > 0) {
    if (verbose) {
      cat("\nModel Comparison:\n")
      cat("  Original (all free) NLL: ", round(fit_all_free$neg_log_lik, 2), "\n")
      cat("  Final (some fixed) NLL:  ", round(fit_final$neg_log_lik, 2), "\n")
      cat("  Original AIC: ", round(fit_all_free$aic, 2), "\n")
      cat("  Final AIC:    ", round(fit_final$aic, 2), "\n")
      cat("  Original BIC: ", round(fit_all_free$bic, 2), "\n")
      cat("  Final BIC:    ", round(fit_final$bic, 2), "\n")
    }
  }
  
  # Return final results
  list(
    fit_all_free = fit_all_free,
    fit_final = fit_final,
    diagnostics = diagnostics,
    params_to_fix = params_to_fix
  )
}

###############################################################################
# Example Usage
###############################################################################

# Generate test data and run full analysis
run_diagnostics_test <- function() {
  # Generate simulation data
  true_params <- list(alpha=0.3, temp=1.5, envy=1.0, guilt=0.5)
  sim_data <- simulate_2games_data_mf_hmm(true_params, game_size=25, playerId=1)
  
  # Run analysis
  results <- analyze_model_identifiability(sim_data, n_recovery=5, verbose=TRUE)
  
  return(results)
}

# Usage:
# test_results <- run_diagnostics_test()