
###### Step 4: Real Data Application ######

source("poisson_ridge.R")

# Load real datasets
wbc <- read.csv("cWBC.csv")
mi_samples <- read.csv("cMISAMPLES.csv")
pid_meta <- read.csv("cPIDMETA.csv")
hct_meta <- read.csv("cHCTMETA.csv")
med <- read.csv("cMED.csv")


str(wbc)
str(mi_samples)
str(pid_meta)
str(hct_meta)
str(med)

names(wbc)
names(mi_samples)
names(pid_meta)
names(hct_meta)
names(med)


###### Preprocessing ######
unique(wbc$wbc)

#### Prepare response variable ####
neutrophils <- subset(
  wbc,
  wbc == "neutrophils"
)

# Days relative to transplantation
neutrophils$day <- neutrophils$anon.date -
  neutrophils$anon.bmt.date

# Distance from Day 14
neutrophils$day_distance <- abs(
  neutrophils$day - 14
)

# Inspect neutrophil data
#dim(neutrophils)
#summary(neutrophils$count)
#summary(neutrophils$day)
#table(neutrophils$exclude_reason, useNA = "ifany")
#sum(is.na(neutrophils$count))
#sum(is.na(neutrophils$day))

# Remove missing counts and restrict the time window
neutrophils_clean <- subset(
  neutrophils,
  !is.na(count) &
    count >= 0 &
    day >= 7 &
    day <= 21 &
    exclude_reason != "cbc_not_closest_to_noon_on_this_day"
)

# Sort observations by distance from Day 14
neutrophils_clean <- neutrophils_clean[
  order(
    neutrophils_clean$pid,
    neutrophils_clean$anon.bmt.date,
    neutrophils_clean$day_distance,
    neutrophils_clean$day
  ),
]

# Select the closest observation for each transplantation
response_data <- neutrophils_clean[
  !duplicated(
    neutrophils_clean[, c("pid", "anon.bmt.date")]
  ),
]


# Keep one transplantation per patient
response_data <- response_data[
  !duplicated(response_data$pid),
]


# Round the response to non-negative integers
response_data$y <- round(response_data$count)


# Check response data
#dim(response_data)
#length(unique(response_data$pid))
#summary(response_data$day)
#summary(response_data$day_distance)
#summary(response_data$y)
#table(response_data$day)
#sum(is.na(response_data$y))
#anyDuplicated(response_data$pid)


#### Prepare Predictors ####

# A. Microbiome Diversity

diversity_data <- mi_samples[
  , c("pid", "anon_date", "inverseSimpson")
]

# Remove missing diversity values
diversity_data <- subset(
  diversity_data,
  !is.na(inverseSimpson) &
    inverseSimpson > 0
)

# Match with response data
diversity_data <- merge(
  diversity_data,
  response_data[, c("pid", "anon.date", "anon.bmt.date")],
  by = "pid"
)

# Calculate time difference
diversity_data$day_distance <-
  diversity_data$anon.date - diversity_data$anon_date

# Keep samples collected within 7 days before response
diversity_data <- subset(
  diversity_data,
  day_distance >= 0 &
    day_distance <= 7 &
    anon_date >= anon.bmt.date
)

# Select the closest microbiome sample
diversity_data <- diversity_data[
  order(
    diversity_data$pid,
    diversity_data$day_distance
  ),
]

diversity_data <- diversity_data[
  !duplicated(diversity_data$pid),
]

# Keep required columns
diversity_data <- diversity_data[
  , c("pid", "inverseSimpson")
]

dim(diversity_data)

summary(diversity_data$inverseSimpson)

anyDuplicated(diversity_data$pid)


# B. Patient Age and Sex

demographic_data <- pid_meta[
  , c("pid", "age", "sex")
]

# Remove identical records
demographic_data <- unique(demographic_data)

# Check remaining duplicate patient IDs
sum(duplicated(demographic_data$pid))

summary(demographic_data$age)

table(demographic_data$sex, useNA = "ifany")

# Inspect conflicting records
head(
  demographic_data[
    duplicated(demographic_data$pid) |
      duplicated(demographic_data$pid, fromLast = TRUE),
  ],
  20
)

# Remove patients with conflicting demographic records

duplicate_pid <- demographic_data$pid[
  duplicated(demographic_data$pid)
]

demographic_data <- subset(
  demographic_data,
  !pid %in% duplicate_pid
)

# Check demographic data

dim(demographic_data)

anyDuplicated(demographic_data$pid)


# C. Graft Type

graft_data <- unique(
  hct_meta[, c("pid", "anon.bmt.date", "hct.source")]
)

graft_data <- merge(
  response_data[, c("pid", "anon.bmt.date")],
  graft_data,
  by = c("pid", "anon.bmt.date")
)

# Check graft data
dim(graft_data)

table(graft_data$hct.source, useNA = "ifany")

anyDuplicated(graft_data$pid)

# D. Drug Administration

head(
  sort(table(med$factor), decreasing = TRUE),
  30
)

# Select drugs

selected_drugs <- c(
  "filgrastim",
  "tacrolimus",
  "hydrocortisone",
  "vancomycin",
  "acyclovir"
)

drug_data <- subset(
  med,
  factor %in% selected_drugs
)


drug_data <- merge(
  drug_data,
  response_data[, c(
    "pid",
    "anon.date",
    "anon.bmt.date"
  )],
  by = "pid"
)

# Keep drug administrations before response measurement

drug_data <- subset(
  drug_data,
  anon.start.date < anon.date &
    anon.stop.date >= anon.bmt.date
)

# Construct drug indicators

drug_indicators <- data.frame(
  pid = response_data$pid
)

for (drug in selected_drugs) {
  
  exposed_pid <- drug_data$pid[
    drug_data$factor == drug
  ]
  
  drug_indicators[[drug]] <- as.integer(
    drug_indicators$pid %in% exposed_pid
  )
}

dim(drug_indicators)

anyDuplicated(drug_indicators$pid)

colSums(drug_indicators[, -1])

colMeans(drug_indicators[, -1])



###### Merge Datasets ######

real_data <- response_data[, c("pid", "y")]

# Merge microbiome diversity
real_data <- merge(
  real_data,
  diversity_data,
  by = "pid"
)

# Merge age and sex
real_data <- merge(
  real_data,
  demographic_data,
  by = "pid"
)

# Merge graft type
real_data <- merge(
  real_data,
  graft_data[, c("pid", "hct.source")],
  by = "pid"
)

# Merge drug administration
real_data <- merge(
  real_data,
  drug_indicators,
  by = "pid"
)

# Remove missing values
real_data <- na.omit(real_data)

dim(real_data)

str(real_data)

summary(real_data)

anyDuplicated(real_data$pid)

colSums(is.na(real_data))

# Remove predictors with zero or near-zero variance
real_data$filgrastim <- NULL
real_data$acyclovir <- NULL


### Standardize continuous predictors ###
real_data$age_std <- as.vector(
  scale(real_data$age)
)

real_data$diversity_std <- as.vector(
  scale(real_data$inverseSimpson)
)


### Convert categorical predictors to factors ###

real_data$sex <- factor(real_data$sex)

real_data$hct.source <- factor(real_data$hct.source)

# Construct design matrix

X_real <- model.matrix(
  ~ diversity_std + age_std + sex + hct.source +
    tacrolimus + hydrocortisone + vancomycin,
  data = real_data
)

# Response vector
y_real <- real_data$y
dim(X_real)
colnames(X_real)
length(y_real)



###### Fit Real Data ######

fit_real <- poisson_ridge(
  X = X_real,
  y = y_real,
  lambda = 0.1,
  epsilon = 1e-3,
  max_iter = 1000,
  rho = 0.5
)

beta_hat_real <- fit_real$beta

# Check convergence
fit_real$converged
fit_real$iterations
fit_real$gradient_norm
fit_real$normalized_gradient_norm


# Check monotonic decrease
all(diff(fit_real$loss_history) <= 0)


plot(
  seq_along(fit_real$loss_history) - 1,
  fit_real$loss_history,
  type = "l",
  col = "blue",
  lwd = 2,
  xlab = "Iteration",
  ylab = "Objective Function",
  main = "Convergence on Real Data"
)

loss_gap_real <- fit_real$loss_history -
  tail(fit_real$loss_history, 1)

plot(
  seq_along(loss_gap_real) - 1,
  pmax(loss_gap_real, 1e-12),
  type = "l",
  log = "y",
  col = "blue",
  lwd = 2,
  xlab = "Iteration",
  ylab = "Objective Gap (log scale)",
  main = "Log-scale Convergence on Real Data"
)

# Estimated coefficients

real_results <- data.frame(
  Parameter = colnames(X_real),
  Estimate = beta_hat_real
)

real_results


##### Poisson Assumption #####
# Predicted conditional means

mu_hat_real <- as.vector(
  exp(X_real %*% beta_hat_real)
)

# Pearson dispersion statistic

dispersion <- sum(
  (y_real - mu_hat_real)^2 / mu_hat_real
) / (nrow(X_real) - ncol(X_real))

dispersion



###### Save Real Data Results ######

# Create output folder
dir.create("results", showWarnings = FALSE)

# Save estimated coefficients
write.csv(
  real_results,
  "results/real_data_coefficients.csv",
  row.names = FALSE
)

# Save convergence results
convergence_results_real <- data.frame(
  n = nrow(X_real),
  p = ncol(X_real) - 1,
  Lambda = 0.1,
  Converged = fit_real$converged,
  Iterations = fit_real$iterations,
  Gradient_Norm = fit_real$gradient_norm,
  Normalized_Gradient_Norm = fit_real$normalized_gradient_norm,
  Monotonic_Decrease = all(diff(fit_real$loss_history) <= 0),
  Final_Objective = tail(fit_real$loss_history, 1),
  Pearson_Dispersion = dispersion
)

write.csv(
  convergence_results_real,
  "results/real_data_summary.csv",
  row.names = FALSE
)

# Save objective function history
loss_results_real <- data.frame(
  Iteration = seq_along(fit_real$loss_history) - 1,
  Objective = fit_real$loss_history,
  Objective_Gap = loss_gap_real
)

write.csv(
  loss_results_real,
  "results/real_data_loss_history.csv",
  row.names = FALSE
)

# Save convergence plot

png(
  "results/real_data_convergence.png",
  width = 1800,
  height = 1200,
  res = 200
)

plot(
  seq_along(fit_real$loss_history) - 1,
  fit_real$loss_history,
  type = "l",
  col = "blue",
  lwd = 2,
  xlab = "Iteration",
  ylab = "Objective Function",
  main = "Convergence on Real Data"
)

dev.off()

# Save log-scale convergence plot

idx <- which(loss_gap_real > 0)

png(
  "results/real_data_log_convergence.png",
  width = 1800,
  height = 1200,
  res = 200
)

plot(
  idx - 1,
  loss_gap_real[idx],
  type = "l",
  log = "y",
  col = "blue",
  lwd = 2,
  xlab = "Iteration",
  ylab = "Objective Gap (log scale)",
  main = "Log-scale Convergence on Real Data"
)

dev.off()

# Save fitted model
saveRDS(
  fit_real,
  "results/real_data_fit.rds"
)

