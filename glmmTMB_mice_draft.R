
library(glmmTMB)
library(mice)
library(miceadds)
library(lme4)



mice.impute.2l.glmmTMB <- function(y,
                                   ry,
                                   x,
                                   type,
                                   wy = NULL,
                                   formula = NULL,
                                   seed = NULL,
                                   ...) {
  
  # Check package
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("Please install glmmTMB first: install.packages('glmmTMB')")
  }
  
  # I keep the seed thing here because it helped me test the function in the first version.
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # wy are rows should be imputed.
  if (is.null(wy)) {
    wy <- !ry
  }
  
  # x a data frame
  x <- as.data.frame(x)
  
  # type must have names
  if (is.null(names(type))) {
    names(type) <- colnames(x)
  }
  
  # -2: the cluster variable.
  #  1: fixed effect.
  #  2: random slope.
  clust <- names(type[type == -2])
  fixe  <- names(type[type > 0])
  rande <- names(type[type == 2])
  
  if (length(clust) != 1) {
    stop("You must have exactly one cluster variable with type = -2.")
  }
  
  dat <- data.frame(y = y, x)
  
  # I make the cluster a factor because it is a grouping variable
  dat[[clust]] <- factor(dat[[clust]])
  
  # Build formula automatically unless the user gives one
  if (!is.null(formula)) {
    randmodel <- as.formula(formula)
  } else {
    
    fixed_part <- if (length(fixe) == 0) {
      "1"
    } else {
      paste(fixe, collapse = " + ")
    }
    
    random_part <- if (length(rande) == 0) {
      "1"
    } else {
      paste("1 +", paste(rande, collapse = " + "))
    }
    
    randmodel <- as.formula(
      paste("y ~", fixed_part, "+ (", random_part, "|", clust, ")")
    )
  }
  
  # Fit model only on observed rows
  fit <- try(
    glmmTMB::glmmTMB(
      formula = randmodel,
      data = dat[ry, ],
      ...
    ),
    silent = TRUE
  )
  
  # If model fails return old missing values
  if (inherits(fit, "try-error")) {
    warning("glmmTMB did not run. I used a simple fallback imputation.")
    return(rnorm(sum(wy), mean = mean(y[ry], na.rm = TRUE), sd = sd(y[ry], na.rm = TRUE)))
  }
  
  # Predict missing values.
  predicted_values <- predict(
    fit,
    newdata = dat[wy, ],
    type = "response",
    allow.new.levels = TRUE
  )
  
  # Adding random noise
  residual_sd <- sigma(fit)
  
  imputed_values <- predicted_values + rnorm(
    n = sum(wy),
    mean = 0,
    sd = residual_sd
  )
  
  return(as.numeric(imputed_values))
}


# -----------------------------------------

set.seed(2026)

# I now make this a simulation study instead of one single simulation
n_sim <- 10

# I save the results of each simulation repetition here
simulation_results_list <- vector("list", n_sim)

for (sim in 1:n_sim) {
  
  # I use a different seed for each simulation repetition
  set.seed(2026 + sim)
  
  cat("\nSimulation repetition:", sim, "of", n_sim, "\n")

# I made the dataset larger than before but we can still make it even larger
# When we make it larger it mice() takes much more time to run the itrerations
n_schools <- 200
students_per_school <- 70
N <- n_schools * students_per_school

school_id <- factor(rep(1:n_schools, each = students_per_school))

age <- rnorm(N, mean = 16, sd = 1)
study_hours <- rnorm(N, mean = 5, sd = 2)

# Random school effects

school_intercept <- rnorm(n_schools, mean = 0, sd = 4)
school_slope <- rnorm(n_schools, mean = 0, sd = 0.4)

u_school <- school_intercept[school_id]
u_slope <- school_slope[school_id]

# True parameter I care about: effect of study_hours on score
true_beta_study_hours <- 1.8

# True score
score_true <- 50 +
  2.5 * age +
  true_beta_study_hours * study_hours +
  u_school +
  u_slope * study_hours +
  rnorm(N, mean = 0, sd = 5)

complete_data <- data.frame(
  school_id = school_id,
  age = age,
  study_hours = study_hours,
  score = score_true
)



# Running analysis on complete data first
# ---------------------------------------

# Before adding missing values, I first fit the analysis model on the full data
# (I still don't understand the purpose but I believe it can give us a reference
# estimate before imputation)

complete_fit <- glmmTMB(
  score ~ age + study_hours + (1 + study_hours | school_id),
  data = complete_data
)

summary(complete_fit)

complete_estimate <- summary(complete_fit)$coefficients$cond["study_hours", "Estimate"]
complete_estimate



# Creating missing values
# Now I am putting missing values mainly in the predictors and also some in the outcome.
# we have a MCAR situation

age_missing <- age
study_hours_missing <- study_hours
score_missing <- score_true

missing_age_rows <- sample(1:N, size = 0.20 * N)
missing_study_rows <- sample(1:N, size = 0.25 * N)
missing_score_rows <- sample(1:N, size = 0.10 * N)

age_missing[missing_age_rows] <- NA
study_hours_missing[missing_study_rows] <- NA
score_missing[missing_score_rows] <- NA

missing_data <- data.frame(
  school_id = school_id,
  age = age_missing,
  study_hours = study_hours_missing,
  score = score_missing
)

# Making sure the cluster variable is integer (I got some errors in imputation because I didn't do this step)
missing_data$school_id <- as.integer(missing_data$school_id)

# Quick check of how many missing values I created
colSums(is.na(missing_data))


# Now I use the function inside mice() with "method =" argument

# We do 5 for now
m <- 5
maxit <- 5

method_glmmTMB <- c(
  school_id = "",
  age = "2l.glmmTMB",
  study_hours = "2l.glmmTMB",
  score = "2l.glmmTMB"
)

# For lmer, I use the already existing method ine miceadds for continuous vars imputation

method_lmer <- c(
  school_id = "",
  age = "2l.continuous",
  study_hours = "2l.continuous",
  score = "2l.continuous"
)

# I now write the prediction matrix manually because I want to clearly see what each model uses
pred <- matrix(
  0,
  nrow = ncol(missing_data),
  ncol = ncol(missing_data),
  dimnames = list(names(missing_data), names(missing_data))
)

# Imputation model for age
pred["age", "school_id"] <- -2
pred["age", "study_hours"] <- 1
pred["age", "score"] <- 1

# Imputation model for study_hours
pred["study_hours", "school_id"] <- -2
pred["study_hours", "age"] <- 1
pred["study_hours", "score"] <- 1

# Imputation model for score
# I keep study_hours as type 2 because I want a random slope for it
pred["score", "school_id"] <- -2
pred["score", "age"] <- 1
pred["score", "study_hours"] <- 2

pred


# Running mice with my glmmTMB function method
# -----------------------------------

# I am not sure if it's necessary to put a value for the seed argument here 
runtime_mice_glmmTMB <- system.time({
  imp_glmmTMB <- mice(
    missing_data,
    m = m,
    maxit = maxit,
    method = method_glmmTMB,
    predictorMatrix = pred,
    printFlag = TRUE,
    seed = 2026
  )
})

runtime_mice_glmmTMB



# Running mice with the existing lmer method
# -----------------------------------

runtime_mice_lmer <- system.time({
  imp_lmer <- mice(
    missing_data,
    m = m,
    maxit = maxit,
    method = method_lmer,
    predictorMatrix = pred,
    printFlag = TRUE,
    seed = 2026
  )
})

runtime_mice_lmer



# Analyzing the imputed datasets from glmmTMB
# ------------------------------------------------------------------------------

# Now I check bias in the later analysis, not only bias in imputed values as I did in the previous code
# To do that I fit the same analysis model on each completed dataset by mice (which used
# glmmTMB as imputation model)
# The parameter I care about for now is the beta coefficient of study_hours.

# I save the estimate coeff and SE here because I need them later for Rubin's rules 
q_glmmTMB <- rep(NA, m)
se_glmmTMB <- rep(NA, m)

# In this loop I analyze each one of the completed datasets created by mice 
for (k in 1:m) {
  completed_k <- complete(imp_glmmTMB, k)
  
  # transformng it back to factor to fit glmmTMB (analysis model)
  completed_k$school_id <- factor(completed_k$school_id)
  
  fit_k <- glmmTMB(
    score ~ age + study_hours + (1 + study_hours | school_id),
    data = completed_k
  )
  
  q_glmmTMB[k] <- summary(fit_k)$coefficients$cond["study_hours", "Estimate"]
  se_glmmTMB[k] <- summary(fit_k)$coefficients$cond["study_hours", "Std. Error"]
}

# After fitting the model on each imputed dataset, I need to combine the m estimates into one final estimate.

# Rubin's rules (written directly so I can understand what happens)
qbar_glmmTMB <- mean(q_glmmTMB) # average estimate
ubar_glmmTMB <- mean(se_glmmTMB^2) # average within-imputation variance
b_glmmTMB <- var(q_glmmTMB) # average between-imputation variance
total_var_glmmTMB <- ubar_glmmTMB + (1 + 1 / m) * b_glmmTMB # total variance (+ correction)
total_se_glmmTMB <- sqrt(total_var_glmmTMB) # pooled standard error

# 95% confidence interval around the pooled estimate:
lower_glmmTMB <- qbar_glmmTMB - 1.96 * total_se_glmmTMB
upper_glmmTMB <- qbar_glmmTMB + 1.96 * total_se_glmmTMB

# Bias compared with the true value tells us whether the multiple imputation 
# overestimates or underestimates the real parameter used in the simulation
bias_true_glmmTMB <- qbar_glmmTMB - true_beta_study_hours

# now we check whether the 95% confidence interval includes the true value
covered_glmmTMB <- lower_glmmTMB <= true_beta_study_hours & upper_glmmTMB >= true_beta_study_hours


# Analyzing the imputed datasets from miceadds lmer 
# (same procedure as for glmmTMB mice function)
# Basically using glmmTMB as analysis model for completed datasets that were imputed using
# lmer as imputation model (we can discuss if this is the correct way to do it)
# ------------------------------------------------------------------------------

q_lmer <- rep(NA, m)
se_lmer <- rep(NA, m)

for (k in 1:m) {
  completed_k <- complete(imp_lmer, k)
  completed_k$school_id <- factor(completed_k$school_id)
  
  fit_k <- glmmTMB(
    score ~ age + study_hours + (1 + study_hours | school_id),
    data = completed_k
  )
  
  q_lmer[k] <- summary(fit_k)$coefficients$cond["study_hours", "Estimate"]
  se_lmer[k] <- summary(fit_k)$coefficients$cond["study_hours", "Std. Error"]
}

qbar_lmer <- mean(q_lmer)
ubar_lmer <- mean(se_lmer^2)
b_lmer <- var(q_lmer)
total_var_lmer <- ubar_lmer + (1 + 1 / m) * b_lmer
total_se_lmer <- sqrt(total_var_lmer)

lower_lmer <- qbar_lmer - 1.96 * total_se_lmer
upper_lmer <- qbar_lmer + 1.96 * total_se_lmer

bias_true_lmer <- qbar_lmer - true_beta_study_hours
bias_complete_lmer <- qbar_lmer - complete_estimate
covered_lmer <- lower_lmer <= true_beta_study_hours & upper_lmer >= true_beta_study_hours



# Evaluation of imputation performance
# ------------------------------------------------------------------------------


# The parameter of interest is the effect of study_hours, the true value is 1.8
# Combing the evaluation into one single dataframe:
# I save one row per method for this simulation repetition
evaluation_sim <- data.frame(
  simulation = sim,
  
  method = c("2l.glmmTMB", "2l.continuous (lmer)"),
  
  estimate = c(qbar_glmmTMB, qbar_lmer),
  
  true_value = true_beta_study_hours,
  
  bias = c(
    qbar_glmmTMB - true_beta_study_hours,
    qbar_lmer - true_beta_study_hours
  ),
  
  absolute_bias = abs(c(
    qbar_glmmTMB - true_beta_study_hours,
    qbar_lmer - true_beta_study_hours
  )),
  
  lower = c(lower_glmmTMB, lower_lmer),
  upper = c(upper_glmmTMB, upper_lmer),
  
  covered = c(covered_glmmTMB, covered_lmer),
  
  ci_width = c(
    upper_glmmTMB - lower_glmmTMB,
    upper_lmer - lower_lmer
  ),
  
  runtime_seconds = c(
    as.numeric(runtime_mice_glmmTMB["elapsed"]),
    as.numeric(runtime_mice_lmer["elapsed"])
  )
)

# I save the result of this repetition
simulation_results_list[[sim]] <- evaluation_sim

print(evaluation_sim)
  }

# Combining all simulation repetitions
evaluation <- do.call(rbind, simulation_results_list)

evaluation


# Summary across the 10 simulation repetitions
simulation_summary <- aggregate(
  cbind(estimate, bias, absolute_bias, covered, ci_width, runtime_seconds) ~ method,
  data = evaluation,
  FUN = mean
)

simulation_summary

# Graph 1: bias

# Smaller absolute bias means the estimate is closer to the true value

barplot(
  evaluation$absolute_bias,
  names.arg = evaluation$method,
  main = "Absolute bias",
  ylab = "Absolute bias"
)


# Graph 2: Confidence interval coverage

# With only a one simulation this is just a rough check not a stable 95% coverage result (to be discussed)

plot(
  1:2,
  evaluation$estimate,
  ylim = range(c(evaluation$lower, evaluation$upper, true_beta_study_hours)),
  xaxt = "n",
  xlab = "Method",
  ylab = "Estimate for study_hours",
  main = "95% CI coverage of the true value"
)

axis(1, at = 1:2, labels = evaluation$method)

arrows(
  x0 = 1:2,
  y0 = evaluation$lower,
  x1 = 1:2,
  y1 = evaluation$upper,
  angle = 90,
  code = 3,
  length = 0.05
)

points(
  1:2,
  evaluation$estimate,
  pch = 19
)

abline(h = true_beta_study_hours, lty = 2, lwd = 2)


# Graph 3: Runtime

barplot(
  evaluation$runtime_seconds,
  names.arg = evaluation$method,
  main = "Runtime comparison",
  ylab = "Seconds"
)


# End of code
