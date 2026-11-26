# Ridge-Penalized Poisson Regression

STAT 607 Final Project

## 1. Introduction

In this project, I implemented ridge-penalized Poisson regression using gradient descent. I used backtracking line search to choose the step size at each iteration.

After implementing the algorithm in R, I tested it on simulated data to check its convergence and how well it recovered the true regression coefficients. I also conducted two scaling experiments and applied the model to a real dataset.

## 2. Files

- `poisson_ridge.R`: Implementation of the algorithm.
- `simulation.R`: Testing the algorithm on simulated data.
- `scaling_experiment.R`: Testing how runtime changes with sample size and number of predictors.
- `real_data.R`: Data preprocessing and model fitting on real data.
- `Model_Testing_Result_Figures/`: Figures from simulation and scaling experiments.
- `Real_Data_Results/`: Results and figures from the real data application.

## 3. How to Use the Algorithm

First, download the files and set the working directory to the folder containing the R scripts.

Load the function:

```r
source("poisson_ridge.R")
```

### Example: Simulated Data

Generate 1,000 observations with 3 predictors.

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

Fit the model using the algorithm:

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

Here, `lambda` is the ridge penalty parameter, `epsilon` is the convergence tolerance, and `rho` is used to reduce the step size during backtracking.

## 4. Example Results

### Convergence

Check whether the algorithm converged:

```r
fit$converged
fit$iterations
fit$normalized_gradient_norm
```

Plot the objective function at each iteration:

```r
plot(
  seq_along(fit$loss_history) - 1,
  fit$loss_history,
  type = "l",
  xlab = "Iteration",
  ylab = "Objective Function"
)
```

Check whether the objective function decreases:

```r
all(diff(fit$loss_history) <= 0)
```

### Ground Truth Recovery

Compare the estimated coefficients with the true coefficients:

```r
data.frame(
  True = beta_true,
  Estimated = fit$beta
)
```

An example from my simulation results:

| Parameter | True | Estimated |
|---|---:|---:|
| beta_0 | 1.000 | 1.000 |
| beta_1 | 0.500 | 0.485 |
| beta_2 | -0.300 | -0.318 |
| beta_3 | 0.800 | 0.788 |

More results can be found in [Model_Testing_Result_Figures](Model_Testing_Result_Figures/).

## 5. Scaling Experiments

I conducted two experiments:

- Change the sample size while keeping the number of predictors fixed at 10.
- Change the number of predictors while keeping the sample size fixed at 1,000.

I compared both total runtime and runtime per iteration because the number of iterations needed for convergence was different across experiments.

To run these experiments:

```r
source("scaling_experiment.R")
```

## 6. Real Data Application

The real dataset is from the following paper:

[The gut microbiota is associated with immune cell dynamics in humans](https://www.nature.com/articles/s41586-020-2971-8)

I used microbiome diversity, patient age, sex, graft type, and drug administration as predictors. The response variable was the neutrophil count closest to Day 14 after transplantation.

The following datasets are needed:

- `cWBC.csv`
- `cMISAMPLES.csv`
- `cPIDMETA.csv`
- `cHCTMETA.csv`
- `cMED.csv`

Download these datasets from the original study and place them in the working directory.

Then run:

```r
source("real_data.R")
```

The model converged on the real data. However, the Pearson dispersion statistic was approximately 3.12, suggesting that the data had overdispersion. This is a limitation of using Poisson regression for this dataset.

The results are saved in [Real_Data_Results](Real_Data_Results/).
