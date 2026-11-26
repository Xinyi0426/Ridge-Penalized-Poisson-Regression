
####### Step 4: Simulation ########

source("poisson_ridge.R")
dir.create("figures", showWarnings = FALSE)

set.seed(607)

n <- 1000
p <- 3

beta_true <- c(1, 0.5, -0.3, 0.8)

# Generate predictors
X <- matrix(rnorm(n * p), nrow = n, ncol = p)

# Add intercept
X <- cbind(1, X)

# True conditional means
mu_true <- as.vector(exp(X %*% beta_true))

# Poisson response
y <- rpois(n, lambda = mu_true)


# Fit the model (based on previous made model)
fit <- poisson_ridge(
  X = X,
  y = y,
  lambda = 0.1,
  epsilon = 1e-3,
  max_iter = 1000,
  rho = 0.5
)

beta_hat <- fit$beta


### Chekck Conv. ###
fit$converged

fit$iterations

fit$gradient_norm


# Conv. plot

png("figures/convergence.png", width = 1000, height = 700)
plot(
  seq_along(fit$loss_history) - 1,
  fit$loss_history,
  type = "l",
  col = "blue",
  lwd = 2,
  xlab = "Iteration",
  ylab = "Objective Function",
  main = "Convergence of Gradient Descent"
)

dev.off()


# Monotonic Decrease
all(diff(fit$loss_history) <= 0)


# Objective gap on log scale

# log_{10}[J(beta_k)-J(beta_final)]
loss_gap <- fit$loss_history - tail(fit$loss_history, 1)

png("figures/log_scale_convergence.png", width = 1000, height = 700)
plot(
  seq_along(loss_gap) - 1,
  pmax(loss_gap, 1e-12),
  type = "l",
  log = "y",
  col = "blue",
  lwd = 2,
  xlab = "Iteration",
  ylab = "Objective Gap (log scale)",
  main = "Log-scale Convergence"
)

dev.off()


#### Assess Ground Truth Recovery ###
beta_hat <- fit$beta

comparison <- data.frame(
  True = beta_true,
  Estimated = beta_hat,
  Error = beta_hat - beta_true
)

comparison


# RMSE
rmse <- sqrt(mean((beta_hat - beta_true)^2))

rmse


# Ground Truth Recovery Plot

png("figures/ground_truth_recovery.png", width = 1000, height = 700)
plot(
  beta_true,
  beta_hat,
  pch = 19,
  col = "blue",
  xlab = "True Coefficients",
  ylab = "Estimated Coefficients",
  main = "Ground Truth Recovery",
  xlim = range(c(beta_true, beta_hat)),
  ylim = range(c(beta_true, beta_hat)),
  asp = 1
)

abline(a = 0, b = 1, col = "red", lty = 2)

dev.off()


# Save estimated coefficients
write.csv(
  comparison,
  "Model_Testing_Result_Figures/Simulation/estimated_coefficients.csv",
  row.names = FALSE
)

# Save RMSE
write.csv(
  data.frame(RMSE = rmse),
  "Model_Testing_Result_Figures/Simulation/rmse.csv",
  row.names = FALSE
)
