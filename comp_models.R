## ----setup, include=FALSE---------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(optimx)

#' Simple Q-Learning model for trust game
#' @param params Vector of parameters: c(alpha, temp)
#' @param data Dataframe with columns: investment, return
#' @return Negative log-likelihood
simple_qlearning <- function(params, data) {
    # Extract parameters
    alpha <- params[1]      # Learning rate
    temp <- params[2]       # Temperature
    
    # Initialize variables
    n_trials <- nrow(data)
    log_lik <- 0
    
    # Define return proportions to learn values for (0 to 1 in 0.167 steps)
    return_props <- seq(0, 1, by = 0.167)
    Q_values <- rep(0, length(return_props))
    
    for(t in 1:n_trials) {
        investment <- data$investment[t]
        actual_return <- data$return[t]
        actual_prop <- actual_return / (3 * investment)
        
        # Get probabilities for each return proportion
        utilities <- Q_values
        utilities_scaled <- utilities - max(utilities)  # For numerical stability
        probabilities <- exp(utilities_scaled / temp)
        probabilities <- probabilities / sum(probabilities)
        
        # Find closest return proportion and update log likelihood
        closest_idx <- which.min(abs(return_props - actual_prop))
        log_lik <- log_lik + log(probabilities[closest_idx])
        
        # Update Q-values if not last trial
        if(t < n_trials) {
            reward <- actual_return  # Simple monetary reward
            prediction_error <- reward - Q_values[closest_idx]
            Q_values[closest_idx] <- Q_values[closest_idx] + alpha * prediction_error
        }
    }
    
    return(-log_lik)  # Return negative log-likelihood for minimization
}

#' Fit Q-learning model to single participant's data
#' @param participant_data Dataframe for one participant
#' @return List with fitted parameters and model metrics
fit_qlearning_single <- function(participant_data) {
    # Initial parameter values
    init_params <- c(alpha = 0.5, temp = 1.0)
    
    # Parameter bounds
    lower <- c(0.01, 0.01)  # alpha, temp
    upper <- c(0.99, 10.0)  # alpha, temp
    
    # Fit model
    fit <- try(optimx(init_params,
                     fn = simple_qlearning,
                     data = participant_data,
                     method = "L-BFGS-B",
                     lower = lower,
                     upper = upper))
    
    if(!inherits(fit, "try-error")) {
        # Get best parameters
        best_fit <- fit[which.min(fit$value), ]
        
        # Calculate AIC and BIC
        n_params <- length(init_params)
        n_obs <- nrow(participant_data)
        neg_log_lik <- best_fit$value
        aic <- 2 * neg_log_lik + 2 * n_params
        bic <- 2 * neg_log_lik + n_params * log(n_obs)
        
        return(list(
            alpha = best_fit$alpha,
            temperature = best_fit$temp,
            neg_log_lik = neg_log_lik,
            aic = aic,
            bic = bic,
            convergence = best_fit$convcode
        ))
    } else {
        return(list(
            alpha = NA,
            temperature = NA,
            neg_log_lik = NA,
            aic = NA,
            bic = NA,
            convergence = NA
        ))
    }
}

#' Fit Q-learning model to all participants
#' @param data Full dataset
#' @return Dataframe with results for each participant
fit_qlearning_all <- function(data) {
    # Get unique participants
    participants <- unique(data$playerId)
    
    # Initialize results dataframe
    results <- data.frame(
        playerId = character(),
        game = character(),
        d_level = character(),
        alpha = numeric(),
        temperature = numeric(),
        neg_log_lik = numeric(),
        aic = numeric(),
        bic = numeric(),
        convergence = numeric(),
        stringsAsFactors = FALSE
    )
    
    # Fit model for each participant and game
    for(pid in participants) {
        # Get participant's games
        participant_games <- data %>%
            filter(playerId == pid) %>%
            group_by(gameNum.f) %>%
            group_split()
        
        # Get D-factor level for this participant
        d_level <- unique(data$d_level[data$playerId == pid])
        
        # Fit model for each game
        for(game_data in participant_games) {
            game_name <- unique(game_data$gameNum.f)
            
            # Prepare game data
            game_data <- game_data %>%
                select(investment, return)
            
            fit <- fit_qlearning_single(game_data)
            
            results <- rbind(results, data.frame(
                playerId = pid,
                game = game_name,
                d_level = d_level,
                alpha = fit$alpha,
                temperature = fit$temperature,
                neg_log_lik = fit$neg_log_lik,
                aic = fit$aic,
                bic = fit$bic,
                convergence = fit$convergence
            ))
        }
    }
    
    return(results)
}




## ---------------------------------------------------------------------------------------------------------------------


# Fit model to all participants
qlearning_results <- fit_qlearning_all(data)

# Print summary statistics
print("Summary of parameter estimates:")
print(summary(qlearning_results[c("alpha", "temperature")]))

print("\nModel fit metrics:")
print(summary(qlearning_results[c("neg_log_lik", "aic", "bic")]))

# Save results
write.csv(qlearning_results, "qlearning_results.csv", row.names = FALSE)

# Create visualizations
library(ggplot2)

# Learning rate by D-factor and game
p1 <- ggplot(qlearning_results, aes(x = game, y = alpha, fill = d_level)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = "Learning Rates by Game and D-Factor Level",
         x = "Game",
         y = "Learning Rate (α)",
         fill = "D-Factor Level")

# Temperature by D-factor and game
p2 <- ggplot(qlearning_results, aes(x = game, y = temperature, fill = d_level)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = "Temperature by Game and D-Factor Level",
         x = "Game",
         y = "Temperature",
         fill = "D-Factor Level")

# Model fit comparison
p3 <- ggplot(qlearning_results, aes(x = game, y = aic, fill = d_level)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = "Model Fit (AIC) by Game and D-Factor Level",
         x = "Game",
         y = "AIC",
         fill = "D-Factor Level")

print(p1)
print(p2)
print(p3)


## ---------------------------------------------------------------------------------------------------------------------
# Load required libraries
library(tidyverse)
library(optimx)

#########################################
# Constants and Utility Functions
#########################################

# Define the three investment states in the Hidden Markov Model
STATE_INVESTMENTS <- c(4, 11, 17)  # Fixed investment amounts characterizing each state
N_STATES <- 3

#' Calculate payoffs for both players in the trust game
calculate_payoffs <- function(investment, return_amount) {
    trustee_payoff <- 3 * investment - return_amount  # Trustee keeps tripled amount minus return
    investor_payoff <- return_amount                  # Investor receives returned amount
    return(list(trustee = trustee_payoff, investor = investor_payoff))
}

#' Calculate Fehr-Schmidt utility incorporating inequality aversion
calculate_fs_utility <- function(own_payoff, other_payoff, envy, guilt) {
    # Calculate disadvantageous inequality (when other gets more than me)
    disadv_inequity <- max(other_payoff - own_payoff, 0)
    # Calculate advantageous inequality (when I get more than other)
    adv_inequity <- max(own_payoff - other_payoff, 0)
    
    # Fehr-Schmidt utility function
    utility <- own_payoff - 
              envy * disadv_inequity -  # Utility reduction from disadvantageous inequality
              guilt * adv_inequity        # Utility reduction from advantageous inequality
    
    return(utility)
}

#' Calculate transition probabilities between investor states (how truste believes investor transitions as a function of net_return they get)
#' @param current_state Current investor state (1, 2, or 3)
#' @param sensitivity Parameter determining how strongly return rate affects transitions
#' @return Vector of transition probabilities to each state
calculate_transition_probs <- function(current_state, investment, return_amount, sensitivity) {
    probs <- rep(0, N_STATES)  
    
    # Calculate net return for investor
    net_return <- return_amount - investment
    # Could normalize by investment to make comparable across different investment amounts
    # net_return_norm <- net_return/investment
    
    if(current_state == 1) {
        # In lowest state: can only stay or move up
        p_up <- exp(sensitivity * net_return)  # Higher net returns increase probability of moving up
        p_stay <- 1                           
        norm_const <- p_up + p_stay           
        probs[1:2] <- c(p_stay, p_up) / norm_const     
        
    } else if(current_state == 3) {
        # In highest state: can only stay or move down
        p_down <- exp(-sensitivity * net_return)  # Negative net returns increase probability of moving down
        p_stay <- 1                               
        norm_const <- p_down + p_stay            
        probs[2:3] <- c(p_down, p_stay) / norm_const      
        
    } else {
        # In middle state: can move either direction
        p_up <- exp(sensitivity * net_return)    
        p_down <- exp(-sensitivity * net_return)  
        p_stay <- 1                               
        norm_const <- p_up + p_stay + p_down     
        probs[1:3] <- c(p_down, p_stay, p_up) / norm_const  
    }
    
    return(probs)
}


## ---------------------------------------------------------------------------------------------------------------------
# function for literature priors for hierarchical approach
get_literature_priors <- function() {
    list(
        means = c(
            envy = 2.0,        # Typical Fehr-Schmidt envy
            guilt = 0.5,       # Typical Fehr-Schmidt guilt
            temp = 5.0,        # Common in trust games
            sensitivity = 0.5,  # Conservative estimate
            alpha_Q = 0.3      # Typical learning rate
        ),
        sds = c(
            envy = 1.0,
            guilt = 0.3,
            temp = 2.0,
            sensitivity = 0.3,
            alpha_Q = 0.2
        )
    )
}



## ---------------------------------------------------------------------------------------------------------------------
#########################################
# Core Model Function
#########################################

#' Main RL model function for trustee decisions
#' @param params Vector of parameters to be fitted:
#'        params[1]: envy (disadvantageous inequity aversion)
#'        params[2]: guilt (advantageous inequity aversion)
#'        params[3]: temp (temperature parameter for choice)
#'        params[4]: sensitivity (sensitivity to return rates for state transitions)
#'        params[5]: alpha_Q (learning rate for Q-values)
#' @param data Dataframe containing trial data
#' @return Negative log-likelihood of observed decisions
trustee_decision_model <- function(params, data, use_priors = FALSE) {
    # Extract parameters
    envy <- params[1]    
    guilt <- params[2]   
    temp <- params[3]      
    sensitivity <- params[4]  
    alpha_Q <- params[5]     
    
    n_trials <- nrow(data)
    log_lik <- 0  # Initialize log likelihood
    
    # Infer initial state from first investment
    current_state <- which.min(abs(STATE_INVESTMENTS - data$investment[1]))
    
    # Initialize Q-values for each state
    Q_values <- rep(0, N_STATES)
    
    for(t in 1:n_trials) {
        investment <- data$investment[t]
        actual_return <- data$return[t]
        
        # Generate possible return amounts
        possible_returns <- seq(0, 3*investment, by = 1)
        
        # Calculate utilities for each possible return
        utilities <- sapply(possible_returns, function(r) {
            payoffs <- calculate_payoffs(investment, r)
            immediate_utility <- calculate_fs_utility(payoffs$trustee, 
                                                   payoffs$investor,
                                                   envy, guilt)
            
            trans_probs <- calculate_transition_probs(current_state, 
                                                    investment,
                                                    r, 
                                                    sensitivity)
            
            future_value <- sum(trans_probs * Q_values)
            
            return(immediate_utility + future_value)
        })
        
        # Convert utilities to probabilities using softmax
        utilities_scaled <- utilities - max(utilities)
        probabilities <- exp(utilities_scaled / temp)
        probabilities <- probabilities / sum(probabilities)
        
        # Ensure no zero probabilities
        probabilities <- pmax(probabilities, 1e-10)
        probabilities <- probabilities / sum(probabilities)
        
        # Update log likelihood
        log_lik <- log_lik + log(probabilities[actual_return + 1])
        
        # Update Q-values if not last trial
        if(t < n_trials) {
            payoffs <- calculate_payoffs(investment, actual_return)
            reward <- calculate_fs_utility(payoffs$trustee, 
                                        payoffs$investor,
                                        envy, guilt)
            
            next_investment <- data$investment[t + 1]
            next_state <- which.min(abs(STATE_INVESTMENTS - next_investment))
            
            prediction_error <- reward + Q_values[next_state] - Q_values[current_state]
            Q_values[current_state] <- Q_values[current_state] + alpha_Q * prediction_error

            
            current_state <- next_state
        }
    }
    
    # Add prior terms if using hierarchical approach
    if(use_priors) {
        priors <- get_literature_priors()
        prior_terms <- sum(dnorm(
            params,
            mean = priors$means,
            sd = priors$sds,
            log = TRUE
        ))
        return(-(log_lik + prior_terms))  # Return negative log posterior
    }
    
    # Return negative log-likelihood for minimization
    if(!is.finite(log_lik)) {
        print("Warning: Non-finite log likelihood")
        return(1e10)
    }
  
    return(-log_lik)
}


## ---------------------------------------------------------------------------------------------------------------------
fit_participant <- function(participant_data,hierarchical = FALSE) {
    print("Starting fit_participant")
  
    # Set bounds and initial values based on approach
    if(hierarchical) {
        priors <- get_literature_priors()
        init_params <- priors$means
        lower_bounds <- c(
            envy = 0.01,
            guilt = 0.01,
            temp = 0.01,
            sensitivity = 0.01,
            alpha_Q = 0.01
        )
        upper_bounds <- c(
            envy = 6.0,
            guilt = 2.0,
            temp = 15.0,
            sensitivity = 2.0,
            alpha_Q = 0.99
        )
        parscale <- priors$sds
    } else {
      # Initial parameter values - using more conservative values
      init_params <- c(
          envy = 1,       
          guilt = 1,        
          temp = 1,        
          sensitivity = 0.1,  
          alpha_Q = 0.1      
      )
      # Parameter bounds
      #Standard ranges in literature are typically 0-4 for envy and 0-1 for guilt
      lower_bounds <- c(0.01, 0.01, 0.01, 0.01, 0.01)
      upper_bounds <- c(6, 2, 15, 2, 0.99)
    }

    # Test function evaluation at initial point
    test_val <- trustee_decision_model(init_params, participant_data)
    print(paste("Initial function value:", test_val))
    
    if(is.finite(test_val)) {
        # Try optimization
        fit <- try(optimx(init_params,
                         fn = trustee_decision_model,
                         data = participant_data,
                         use_priors = hierarchical,  # New parameter
                         method = "L-BFGS-B",
                         lower = lower_bounds,
                         upper = upper_bounds,
                         control = list(
                             maxit = 1000,
                             parscale = c(1, 1, 1, 1, 1)
                         )))
        
        if(!inherits(fit, "try-error")) {
            best_fit <- fit[which.min(fit$value), , drop=FALSE]
            
            n_params <- length(init_params)
            n_obs <- nrow(participant_data)
            aic <- 2 * best_fit$value + 2 * n_params
            bic <- 2 * best_fit$value + n_params * log(n_obs)
            
            return(list(
                parameters = c(
                    envy = best_fit$envy,
                    guilt = best_fit$guilt,
                    temp = best_fit$temp,
                    sensitivity = best_fit$sensitivity,
                    alpha_Q = best_fit$alpha_Q
                ),
                log_lik = -best_fit$value,
                aic = aic,
                bic = bic,
                convergence = best_fit$convcode
            ))
        }
    }
    
    return(list(
        parameters = c(envy = NA, guilt = NA, temp = NA, sensitivity = NA, alpha_Q = NA),
        log_lik = NA,
        aic = NA,
        bic = NA,
        convergence = NA
    ))
}


## ---------------------------------------------------------------------------------------------------------------------
# # Test with one participant first
# test_single <- function(hierarchical = FALSE) {
#     first_participant <- unique(final_data$playerId)[1]
#     test_data_single <- final_data %>%
#         dplyr::filter(playerId == first_participant)
#     
#     result <- fit_participant(test_data_single,hierarchical)
#     print("Single participant result:")
#     print(result)
#     
#     return(result)
# }
# 
# # Run single test
# single_result <- test_single()
# single_result2 <- test_single(hierarchical = TRUE)




## ---------------------------------------------------------------------------------------------------------------------
fit_all_participants <- function(data, hierarchical = FALSE) {
    # Split data by participant and game condition
    participants <- split(data, list(data$playerId, data$gameNum.f))
    
    # Add error checking
    if(length(participants) == 0) {
        stop("No participants found in data")
    }

    # Export needed functions and objects to the cluster again to ensure availability
    clusterExport(cl, c("calculate_payoffs", 
                       "calculate_fs_utility",
                       "calculate_transition_probs", 
                       "trustee_decision_model",
                       "fit_participant", 
                       "STATE_INVESTMENTS", 
                       "N_STATES"))
    
    # Load required packages on each node
    clusterEvalQ(cl, {
        library(optimx)
    })
    
    # Use tryCatch to handle errors in parallel processing
    results <- tryCatch({
        # Process participants sequentially if small number, otherwise parallel
        if(length(participants) <= 4) {
            results <- lapply(participants, function(p) {
                tryCatch(fit_participant(p, hierarchical = hierarchical),
                        error = function(e) {
                            warning(paste("Error fitting participant:", e$message))
                            return(NULL)
                        })
            })
        } else {
            results <- parLapplyLB(cl, participants, function(p) {
                tryCatch(fit_participant(p, hierarchical = hierarchical),
                        error = function(e) {
                            warning(paste("Error fitting participant:", e$message))
                            return(NULL)
                        })
            })
        }
        results
    }, error = function(e) {
        warning(paste("Error in parallel processing:", e$message))
        return(NULL)
    })
    
    # Convert results to data frame with better error handling
    results_df <- do.call(rbind, lapply(seq_along(participants), function(i) {
        name <- names(participants)[i]
        fit <- results[[i]]
        
        if(is.null(fit) || all(is.na(fit$parameters))) {
            return(data.frame(
                participant = strsplit(name, "\\.")[[1]][1],
                game = strsplit(name, "\\.")[[1]][2],
                envy = NA, guilt = NA,
                temp = NA, sensitivity = NA, alpha_Q = NA,
                log_lik = NA, 
                aic = NA,           # Added
                bic = NA,           # Added
                convergence = NA,
                error = "Fitting failed",
                approach = ifelse(hierarchical, "hierarchical", "standard")  # Track approach
            ))
        }
        
        data.frame(
            participant = strsplit(name, "\\.")[[1]][1],
            game = strsplit(name, "\\.")[[1]][2],
            envy = fit$parameters["envy"],
            guilt = fit$parameters["guilt"],
            temp = fit$parameters["temp"],
            sensitivity = fit$parameters["sensitivity"],
            alpha_Q = fit$parameters["alpha_Q"],
            log_lik = fit$log_lik,
            aic = fit$aic,           # Added
            bic = fit$bic,           # Added
            convergence = fit$convergence,
            error = NA,
            approach = ifelse(hierarchical, "hierarchical", "standard")  # Track approach
        )
    }))
    
    # Add summary statistics for the group
    attr(results_df, "summary") <- list(
        mean_aic = mean(results_df$aic, na.rm = TRUE),
        mean_bic = mean(results_df$bic, na.rm = TRUE),
        mean_log_lik = mean(results_df$log_lik, na.rm = TRUE),
        param_means = colMeans(results_df[, c("envy", "guilt", "temp",
                                            "sensitivity", "alpha_Q")],
                             na.rm = TRUE),
        param_sds = apply(results_df[, c("envy", "guilt", "temp",
                                       "sensitivity", "alpha_Q")],
                         2, sd, na.rm = TRUE),
        convergence_rate = mean(!is.na(results_df$convergence)),
        n_participants = nrow(results_df),
        approach = ifelse(hierarchical, "hierarchical", "standard")
    )

    return(results_df)
}


