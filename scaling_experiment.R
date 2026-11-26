
###### Step 5: Expected Scaling Behavior #######

source("poisson_ridge.R")
dir.create("figures", showWarnings = FALSE)

# Simulation Function
set.seed(607)

generate_data <- function(n, p) {
  
  X <- cbind(1, matrix(rnorm(n * p), nrow = n))
  
  beta_true <- c(1, rep(0.8 / sqrt(p), p))
  
  mu_true <- as.vector(exp(X %*% beta_true))
  
  y <- rpois(n, lambda = mu_true)
  
  return(list(
    X = X,
    y = y,
    beta_true = beta_true
  ))
}


# Runtime Measurement Function

measure_runtime <- function(n, p, lambda = 1,
                            reps = 5) {
  
  runtime <- numeric(reps)
  iterations <- numeric(reps)
  
  for (i in seq_len(reps)) {
    
    data <- generate_data(n, p)
    
    runtime[i] <- system.time({
      
      fit <- poisson_ridge(
        X = data$X,
        y = data$y,
        lambda = lambda,
        epsilon = 1e-6
      )
      
    })["elapsed"]
    
    iterations[i] <- length(fit$loss_history) - 1
  }
  
  return(data.frame(
    n = n,
    p = p,
    Runtime = median(runtime),
    Iterations = median(iterations)
  ))
}


#Experiment 1: Scaling with Sample Size n

n_values <- c(1000, 2000, 5000,
              10000, 20000, 50000)

p_fixed <- 10

results_n <- do.call(rbind, lapply(n_values, function(n) {
  
  measure_runtime(
    n = n,
    p = p_fixed
  )
  
}))

results_n

# Runtime vs. Sample Size

png("figures/scaling_sample_size.png", width = 1000, height = 700)
plot(
  results_n$n,
  results_n$Runtime,
  type = "b",
  pch = 19,
  col = "blue",
  lwd = 2,
  xlab = "Number of Observations (n)",
  ylab = "Runtime (seconds)",
  main = "Scaling with Sample Size"
)

dev.off()


# Runtime per iteration
results_n$Time_per_iter <- 
  results_n$Runtime / results_n$Iterations

results_n

png("figures/scaling_sample_size_per_iteration.png", width = 1000, height = 700)
plot(
  results_n$n,
  results_n$Time_per_iter,
  type = "b",
  pch = 19,
  col = "blue",
  lwd = 2,
  xlab = "Number of Observations (n)",
  ylab = "Runtime per Iteration (seconds)",
  main = "Per-Iteration Scaling with Sample Size"
)

dev.off()


# Experiment 2: Scaling with p

p_values <- c(5, 10, 20, 50, 100, 200, 500)

n_fixed <- 1000

results_p <- do.call(rbind, lapply(p_values, function(p) {
  
  measure_runtime(
    n = n_fixed,
    p = p
  )
  
}))

results_p

# Runtime vs. Number of Parameters

png("figures/scaling_predictors.png", width = 1000, height = 700)
plot(
  results_p$p,
  results_p$Runtime,
  type = "b",
  pch = 19,
  col = "blue",
  lwd = 2,
  xlab = "Number of Predictors (p)",
  ylab = "Runtime (seconds)",
  main = "Scaling with Number of Parameters"
)

dev.off()


# Runtime per iteration

results_p$Time_per_iter <- 
  results_p$Runtime / results_p$Iterations

results_p

png("figures/scaling_predictors_per_iteration.png", width = 1000, height = 700)
plot(
  results_p$p,
  results_p$Time_per_iter,
  type = "b",
  pch = 19,
  col = "blue",
  lwd = 2,
  xlab = "Number of Predictors (p)",
  ylab = "Runtime per Iteration (seconds)",
  main = "Per-Iteration Scaling with Predictors"
)

dev.off()

# Save numerical results

write.csv(
  results_n,
  "Model_Testing_Result_Figures/scaling_results_n.csv",
  row.names = FALSE
)

write.csv(
  results_p,
  "Model_Testing_Result_Figures/scaling_results_p.csv",
  row.names = FALSE
)

