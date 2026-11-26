
###### Step 3: Gradient Descent Modeling ######
# Objective Function
J <- function(beta_k, X, y, lambda, D) {
  
  mu_k <- as.vector(exp(X %*% beta_k))
  
  nll <- sum(mu_k - y * as.vector(X %*% beta_k))
  
  penalty <- (lambda / 2) * t(beta_k) %*% D %*% beta_k
  
  return(as.numeric(nll + penalty))
}



# Gradient Descent
poisson_ridge <- function(X, y, lambda, 
                          epsilon = 1e-3, 
                          max_iter = 1000,
                          rho = 0.5) {
  
  # Set initial values
  k <- 0
  beta_k <- rep(0, ncol(X))
  
  D <- diag(ncol(X))
  D[1, 1] <- 0
  
  loss_history <- numeric()
  
  while (k < max_iter) {
    
    # 1. Predicted mean
    mu_k <- as.vector(exp(X %*% beta_k))
    
    
    # 2. Objective function
    J_k <- J(beta_k, X, y, lambda, D)
    
    loss_history <- c(loss_history, J_k)
    
    
    # 3. Gradient
    g_k <- as.vector(
      t(X) %*% (mu_k - y) + lambda * D %*% beta_k
    )
    
    if (sqrt(sum(g_k^2)) / nrow(X) < epsilon) {
      break
    }
    
    
    # 4. Descent direction
    d_k <- -g_k
    
    
    # 5. Choose a step size
    alpha_k <- 1
    
    J_new <- J(beta_k + alpha_k * d_k, 
               X, y, lambda, D)
    
    while (!is.finite(J_new) || J_new >= J_k) {
      
      alpha_k <- rho * alpha_k
      
      if (alpha_k < 1e-16) {
        
        print(list(
          iteration = k,
          beta_k = beta_k,
          gradient_norm = sqrt(sum(g_k^2)),
          current_loss = J_k,
          new_loss = J_new,
          alpha_k = alpha_k
        ))
        
        stop("Backtracking failed")
      }
      
      J_new <- J(beta_k + alpha_k * d_k, 
                 X, y, lambda, D)
    }
    
    
    # 6. Update the parameters
    beta_k <- beta_k + alpha_k * d_k
    
    
    # 7. Update iteration number
    k <- k + 1
  }
  
  # Save the final objective value if iteration limit is reached
  if (k == max_iter) {
    loss_history <- c(
      loss_history, J(beta_k, X, y, lambda, D)
    )
  }
  
  # Calculate the final gradient
  mu_k <- as.vector(exp(X %*% beta_k))
  
  g_k <- as.vector(
    t(X) %*% (mu_k - y) + lambda * D %*% beta_k
  )
  
  gradient_norm <- sqrt(sum(g_k^2))
  
  normalized_gradient_norm <- gradient_norm / nrow(X)
  
  return(list(
    beta = beta_k,
    loss_history = loss_history,
    iterations = k,
    gradient_norm = gradient_norm,
    normalized_gradient_norm = normalized_gradient_norm,
    converged = normalized_gradient_norm < epsilon
  ))
}
