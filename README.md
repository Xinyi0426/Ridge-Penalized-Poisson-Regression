# Ridge-Penalized Poisson Regression

STAT 607 Final Project

## Overview

In this project, I implemented ridge-penalized Poisson regression in R using gradient descent with backtracking line search.

I tested the algorithm on simulated data to check its convergence and coefficient recovery. I also ran two scaling experiments and applied the algorithm to a real dataset.

## Files

- `poisson_ridge.R`: Algorithm implementation
- `simulation.R`: Simulation and model testing
- `scaling_experiment.R`: Runtime and scalability experiments
- `real_data.R`: Real data preprocessing and model fitting
- `Model_Testing_Result_Figures/`: Simulation and scaling results
- `Real_Data_Results/`: Real data application results

## How to use

Download the files and set the working directory to the project folder.

Load the algorithm:

```r
source("poisson_ridge.R")
```

Generate some simulated data:

```r
set.seed(607)

n <- 1000
p <- 3
beta_true <- c(1, 0.5, -0.3, 0.8)

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
X <- cbind(1, X)

mu_true <- as.vector(exp(X %*% beta_true))
y <- rpois(n, lambda = mu_true)
```

Fit the model:

```r
fit <- poisson_ridge(
  X = X,
  y = y,
  lambda = 0.1,
  epsilon = 1e-6,
  max_iter = 1000,
  rho = 0.5
)
```

Here, `lambda` controls the ridge penalty, `epsilon` is the convergence tolerance, and `rho` controls the step size reduction during backtracking.

Check convergence:

```r
fit$converged
fit$iterations
fit$normalized_gradient_norm
```

Plot the objective function:

```r
plot(
  seq_along(fit$loss_history) - 1,
  fit$loss_history,
  type = "l",
  xlab = "Iteration",
  ylab = "Objective Function"
)
```

Check whether the objective function is non-increasing:

```r
all(diff(fit$loss_history) <= 0)
```

Compare the estimated coefficients with the true coefficients:

```r
data.frame(
  True = beta_true,
  Estimated = fit$beta
)
```

Below is an example from my simulation results:

| Parameter | True | Estimated |
|---|---:|---:|
| beta_0 | 1.000 | 1.000 |
| beta_1 | 0.500 | 0.485 |
| beta_2 | -0.300 | -0.318 |
| beta_3 | 0.800 | 0.788 |

More simulation results can be found in `Model_Testing_Result_Figures/`.

## Scaling experiments

I tested the scalability of the algorithm in two settings:

- Varying sample size with 10 predictors.
- Varying the number of predictors with 1,000 observations.

Since the number of iterations needed for convergence varied, I compared both total runtime and runtime per iteration.

To reproduce the experiments:

```r
source("scaling_experiment.R")
```

## Real data application

The dataset used in this project is from the following study:

[The gut microbiota is associated with immune cell dynamics in humans](https://www.nature.com/articles/s41586-020-2971-8)

I used microbiome diversity, age, sex, graft type, and drug administration as predictors. The response variable was the neutrophil count closest to Day 14 after transplantation.

The following data files are needed:

- `cWBC.csv`
- `cMISAMPLES.csv`
- `cPIDMETA.csv`
- `cHCTMETA.csv`
- `cMED.csv`

Download the data from the original study and place the files in the working directory.

Run:

```r
source("real_data.R")
```

The algorithm converged on the real dataset. However, the Pearson dispersion statistic was approximately 3.12, indicating overdispersion. This is one limitation of applying Poisson regression to this dataset.

The results and figures are available in `Real_Data_Results/`.
