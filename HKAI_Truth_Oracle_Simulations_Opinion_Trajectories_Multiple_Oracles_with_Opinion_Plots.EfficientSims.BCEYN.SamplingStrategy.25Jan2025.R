# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Simulation parameters
timesteps <- 10000
num_individuals <- 500
truth_value <- 0.8
epsilon_values <- c(0.08)
truth_weights <- c(0.04, 0.32)
oracle_weights <- c(0.04, 0.32) # Updated to avoid zero weight
p_oracle_access <- c(0.0,0.05 ,0.5, 1.0)
p_truth_access <- c(0.0,0.05,0.5, 1.0)
num_oracles <- 2
rand_sample_proportion <- 1 #specifies the proportion of opinions that Oracles 1 or 2 have access to when sampling is set at "random"
oracle_sampling_strategy <- c("truth-biased") # Options: "random", "truth-biased"
num_simulations <- 5
convergence_tolerance <- 1e-4 # Numerical tolerance for convergence
Op_BCEYN<-c(0,1)

# Function to calculate the number of opinion groups -- not used in Rodrigo, 2025
calculate_groups <- function(opinions, epsilon) {
  group_count <- 0
  visited <- rep(FALSE, length(opinions))
  
  for (i in seq_along(opinions)) {
    if (!visited[i]) {
      group <- abs(opinions - opinions[i]) <= epsilon
      visited <- visited | group
      group_count <- group_count + 1
    }
  }
  return(group_count)
}
# Simulation function
run_simulation <- function(Strategy, BCEYN, epsilon, truth_weight, oracle_weight, p_oracle_access, p_truth_access, 
                           num_individuals, timesteps, truth_value, tolerance, store_trajectory = FALSE) {
  opinions <- runif(num_individuals, 0, 1)
  oracle_access <- rbinom(num_individuals, 1, p_oracle_access)
  truth_access <- rbinom(num_individuals, 1, p_truth_access)
  oracle_opinions <- rep(0, num_oracles)
  
  time_to_converge <- NA
  trajectory <- if (store_trajectory) data.frame(timestep = integer(), individual = integer(), opinion = numeric()) else NULL
  oracle_trajectory <- if (store_trajectory) data.frame(timestep = integer(), oracle = integer(), opinion = numeric()) else NULL
  
  stop_simulation <- FALSE
  
  for (t in 1:timesteps) {
    if (stop_simulation) break
    
    # Calculate Oracle opinions
    for (o in 1:num_oracles) {
      if (Strategy == "random") {
        sampled_indices <- sample(1:num_individuals, size = round(rand_sample_proportion * num_individuals))
        oracle_opinions[o] <- mean(opinions[sampled_indices], na.rm = TRUE)
      } else if (Strategy == "truth-biased") {
        if (o == 1) {
          truth_indices <- which(truth_access == 1)
          oracle_opinions[o] <- mean(opinions[truth_indices], na.rm = TRUE)
        } else if (o > 1) {
          non_truth_indices <- which(truth_access == 0)
          oracle_opinions[o] <- mean(opinions[non_truth_indices], na.rm = TRUE)
        }
      } else {
        stop("Invalid oracle_sampling_strategy or incompatible number of Oracles")
      }
    }
    
    # Update individual opinions
    new_opinions <- numeric(num_individuals)
    for (i in 1:num_individuals) {
      oracle_contrib <- 0
      total_oracle_weight <- 0
      
      for (o in 1:num_oracles) {
        
        if (!is.na(oracle_opinions[o]) && oracle_access[i] == 1 && BCEYN==0) {   #THIS ALLOWS INDIVIDUALS TO ACCESS ALL ORACLES, EVEN IF ITS OUTSIDE THE BCE
          oracle_contrib <- oracle_contrib + oracle_weight * oracle_opinions[o]
          total_oracle_weight <- total_oracle_weight + 1
        } else if (!is.na(oracle_opinions[o]) && oracle_access[i] == 1 && BCEYN==1 && abs(oracle_opinions[o] - opinions[i]) <= epsilon) {   #THIS IS AN IMPORTANT LINE -- IT ALLOWS INDIVIDUALS TO ACCESS AN ORACLE ONLY IF THE ORACLE'S OPINIONS ARE WITHIN THE BCE: 
            oracle_contrib <- oracle_contrib + oracle_weight * oracle_opinions[o]
            total_oracle_weight <- total_oracle_weight + 1
        } else if (oracle_access[i] == 0 ) {    
          oracle_contrib <- 0
          total_oracle_weight <- 0
        }
      }
      # Normalize Oracle contributions if total weight is non-zero
      if (total_oracle_weight > 0) {
        oracle_contrib <- oracle_contrib / total_oracle_weight
      } else {
        oracle_contrib <- 0
        total_oracle_weight <- 0
      }
      
      # Truth contribution
      truth_contrib <- if (truth_access[i] == 1) truth_weight * truth_value else 0
      
      # Bounded Confidence Envelope (BCE) contributions
      bce_indices <- which(abs(opinions - opinions[i]) <= epsilon)
      weighted_bce <- sum(opinions[bce_indices]) / length(bce_indices)
      
      # Remaining weight for Truth and BCE
      remaining_weight <- if (oracle_contrib == 0) 1 else 1 - oracle_weight
      
      # Update opinion
      new_opinion <- if (truth_access[i] == 1) oracle_contrib + remaining_weight * (truth_contrib + (1 - truth_weight) * weighted_bce) else oracle_contrib + remaining_weight * weighted_bce
      
      # No clipping
      new_opinions[i] <- new_opinion
    }
    
    if (is.na(time_to_converge) && all(abs(new_opinions - opinions) < tolerance)) {
      time_to_converge <- t
      stop_simulation <- FALSE
    }
    
    if (!is.na(time_to_converge) && t >= time_to_converge + 10) {
      stop_simulation <- TRUE
    }
    
    if (store_trajectory) {
      trajectory <- rbind(trajectory, data.frame(timestep = t, individual = 1:num_individuals, opinion = new_opinions))
      oracle_trajectory <- rbind(oracle_trajectory, data.frame(timestep = t, oracle = 1:num_oracles, opinion = oracle_opinions))
    }
    
    opinions <- new_opinions
  }
  
  final_variance <- var(opinions)
  group_count <- calculate_groups(opinions, epsilon)
  polarization <- sd(opinions)
  truth_alignment <- mean(abs(opinions - truth_value))
  
  return(list(
    variance = final_variance, 
    groups = group_count, 
    convergence_time = time_to_converge,
    polarization = polarization,
    truth_alignment = truth_alignment,
    trajectory = trajectory,
    oracle_trajectory = oracle_trajectory
  ))
}


# Run batch simulations
results <- data.frame()
trajectories <- data.frame()
oracle_trajectories <- data.frame()

for (Strategy in oracle_sampling_strategy){
for (BCEYN in Op_BCEYN){
for (oracle_weight in oracle_weights) {
  for (epsilon in epsilon_values) {
    for (truth_weight in truth_weights) {
      for (p_oracle in p_oracle_access) {
        for (p_truth in p_truth_access) {
          for (sim in 1:num_simulations) {
            store_trajectory <- (num_simulations == 1)
            result <- run_simulation(
              Strategy, BCEYN, epsilon, truth_weight, oracle_weight, p_oracle, p_truth, 
              num_individuals, timesteps, truth_value, convergence_tolerance, store_trajectory
            )
            
            results <- rbind(results, data.frame(
              Strategy=Strategy,
              BCEYN=BCEYN,
              epsilon = epsilon,
              truth_weight = truth_weight,
              oracle_weight = oracle_weight,
              p_oracle_access = p_oracle,
              p_truth_access = p_truth,
              variance = result$variance,
              groups = result$groups,
              convergence_time = result$convergence_time,
              polarization = result$polarization,
              truth_alignment = result$truth_alignment
            ))
            
            if (!is.null(result$trajectory)) {
              trajectories <- rbind(trajectories, result$trajectory)
            }
            
            if (!is.null(result$oracle_trajectory)) {
              oracle_trajectories <- rbind(oracle_trajectories, result$oracle_trajectory)
            }
          }
        }
      }
    }
  }
}
}
}
# Export results
write.csv(results, "simulation_results.csv", row.names = FALSE)
if (num_simulations == 1) {
  write.csv(trajectories, "opinion_trajectories.csv", row.names = FALSE)
  write.csv(oracle_trajectories, "oracle_trajectories.csv", row.names = FALSE)
}

# Define parameter sets for specific trajectory plots -- use this to produce plots that you want
parameter_sets <- list(
  SET1  = list(Strategy="truth-biased", BCEYN=0, epsilon = 0.08, truth_weight = 0.04, oracle_weight = 0.32, p_oracle_access = 1.0, p_truth_access = 0.05),
  SET2  = list(Strategy="truth-biased", BCEYN=0,epsilon = 0.08, truth_weight = 0.32, oracle_weight = 0.32, p_oracle_access = 1.0, p_truth_access = 0.05),
  SET3  = list(Strategy="truth-biased", BCEYN=0, epsilon = 0.08, truth_weight = 0.04, oracle_weight = 0.32, p_oracle_access = 1, p_truth_access = 0.5),
  SET4  = list(Strategy="truth-biased", BCEYN=0, epsilon = 0.08, truth_weight = 0.32, oracle_weight = 0.32, p_oracle_access = 1, p_truth_access = 0.5)
  #SET5  = list(epsilon = 0.02, truth_weight = 0.01, oracle_weight = 0.32, p_oracle_access = 0.0, p_truth_access = 0.0),
  #SET6  = list(epsilon = 0.02, truth_weight = 0.01, oracle_weight = 0.32, p_oracle_access = 0.0, p_truth_access = 1.0),
  #SET7  = list(epsilon = 0.02, truth_weight = 0.01, oracle_weight = 0.32, p_oracle_access = 1.0, p_truth_access = 0.0),
  #SET8  = list(epsilon = 0.02, truth_weight = 0.01, oracle_weight = 0.32, p_oracle_access = 1.0, p_truth_access = 1.0),
  #SET9  = list(epsilon = 0.02, truth_weight = 0.32, oracle_weight = 0.01, p_oracle_access = 0.0, p_truth_access = 0.0),
  #SET10 = list(epsilon = 0.02, truth_weight = 0.32, oracle_weight = 0.01, p_oracle_access = 0.0, p_truth_access = 1.0),
  #SET11 = list(epsilon = 0.02, truth_weight = 0.32, oracle_weight = 0.01, p_oracle_access = 1.0, p_truth_access = 0.0),
  #SET12 = list(epsilon = 0.02, truth_weight = 0.32, oracle_weight = 0.01, p_oracle_access = 1.0, p_truth_access = 1.0),
  #SET13 = list(epsilon = 0.02, truth_weight = 0.32, oracle_weight = 0.32, p_oracle_access = 0.0, p_truth_access = 0.0),
  #SET14 = list(epsilon = 0.02, truth_weight = 0.32, oracle_weight = 0.32, p_oracle_access = 0.0, p_truth_access = 1.0),
  #SET15 = list(epsilon = 0.02, truth_weight = 0.32, oracle_weight = 0.32, p_oracle_access = 1.0, p_truth_access = 0.0),
  #SET16 = list(epsilon = 0.02, truth_weight = 0.32, oracle_weight = 0.32, p_oracle_access = 1.0, p_truth_access = 1.0),
  #SET17 = list(epsilon = 0.16, truth_weight = 0.01, oracle_weight = 0.01, p_oracle_access = 0.0, p_truth_access = 0.0),
  #SET18 = list(epsilon = 0.16, truth_weight = 0.01, oracle_weight = 0.01, p_oracle_access = 0.0, p_truth_access = 1.0),
  #SET19 = list(epsilon = 0.16, truth_weight = 0.01, oracle_weight = 0.01, p_oracle_access = 1.0, p_truth_access = 0.0),
  #SET20 = list(epsilon = 0.16, truth_weight = 0.01, oracle_weight = 0.01, p_oracle_access = 1.0, p_truth_access = 1.0),
  #SET21 = list(epsilon = 0.16, truth_weight = 0.01, oracle_weight = 0.32, p_oracle_access = 0.0, p_truth_access = 0.0),
  #SET22 = list(epsilon = 0.16, truth_weight = 0.01, oracle_weight = 0.32, p_oracle_access = 0.0, p_truth_access = 1.0),
  #SET23 = list(epsilon = 0.16, truth_weight = 0.01, oracle_weight = 0.32, p_oracle_access = 1.0, p_truth_access = 0.0),
  #SET24 = list(epsilon = 0.16, truth_weight = 0.01, oracle_weight = 0.32, p_oracle_access = 1.0, p_truth_access = 1.0),
  #SET25 = list(epsilon = 0.16, truth_weight = 0.32, oracle_weight = 0.01, p_oracle_access = 0.0, p_truth_access = 0.0),
  #SET26 = list(epsilon = 0.16, truth_weight = 0.32, oracle_weight = 0.01, p_oracle_access = 0.0, p_truth_access = 1.0),
  #SET27 = list(epsilon = 0.16, truth_weight = 0.32, oracle_weight = 0.01, p_oracle_access = 1.0, p_truth_access = 0.0),
  #SET28 = list(epsilon = 0.16, truth_weight = 0.32, oracle_weight = 0.01, p_oracle_access = 1.0, p_truth_access = 1.0),
  #SET29 = list(epsilon = 0.16, truth_weight = 0.32, oracle_weight = 0.32, p_oracle_access = 0.0, p_truth_access = 0.0),
  #SET30 = list(epsilon = 0.16, truth_weight = 0.32, oracle_weight = 0.32, p_oracle_access = 0.0, p_truth_access = 1.0),
  #SET31 = list(epsilon = 0.16, truth_weight = 0.32, oracle_weight = 0.32, p_oracle_access = 1.0, p_truth_access = 0.0),
  #SET32 = list(epsilon = 0.16, truth_weight = 0.32, oracle_weight = 0.32, p_oracle_access = 1.0, p_truth_access = 1.0)
)





# Generate plots for specified parameter sets
for (set_name in names(parameter_sets)) {
  params <- parameter_sets[[set_name]]
  
  # Run simulation for the specific parameter set
  result <- run_simulation(
    Strategy=params$Strategy,
    BCEYN=params$BCEYN,
    epsilon = params$epsilon,
    truth_weight = params$truth_weight,
    oracle_weight = params$oracle_weight,
    p_oracle_access = params$p_oracle_access,
    p_truth_access = params$p_truth_access,
    num_individuals = num_individuals,
    timesteps = timesteps,
    truth_value = truth_value,
    tolerance = convergence_tolerance,
    store_trajectory = TRUE
  )
  
  # Extract trajectories
  opinion_traj <- result$trajectory
  oracle_traj <- result$oracle_trajectory
  
  # Create Opinion Trajectories plot
  opinion_plot <- ggplot(opinion_traj, aes(x = timestep, y = opinion, group = individual, color = as.factor(individual))) +
    geom_line() +
    geom_hline(yintercept = truth_value, linetype = "dashed", color = "black", linewidth = 1) +
    labs(
      title = paste("Opinion Trajectories for", set_name),
      x = "Timestep",
      y = "Opinion",
      color = "Individual"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Create Oracle Trajectories plot
  oracle_plot <- ggplot(oracle_traj, aes(x = timestep, y = opinion, group = oracle, color = as.factor(oracle))) +
    geom_line() +
    geom_hline(yintercept = truth_value, linetype = "dashed", color = "black", linewidth = 1) +
    labs(
      title = paste("Oracle Trajectories for", set_name),
      x = "Timestep",
      y = "Oracle Opinion",
      color = "Oracle"
    ) +
    theme_minimal()
  
  # Print plots side by side
  print(opinion_plot)
  print(oracle_plot)
}
