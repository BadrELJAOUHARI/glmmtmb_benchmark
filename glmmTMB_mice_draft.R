library(glmmTMB)
library(mice)
library(miceadds)
library(lme4)
library(pbapply)


# -----------------------------------------

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
      REML = FALSE,
      ...
    ),
    silent = TRUE
  )
  
  # If model fails return old missing values
  if (inherits(fit, "try-error")) {
    warning("glmmTMB did not run. I used a simple fallback imputation.")
    
    fallback_sd <- sd(y[ry], na.rm = TRUE)
    
    if (is.na(fallback_sd) || fallback_sd == 0) {
      fallback_sd <- 1
    }
    
    return(
      rnorm(
        sum(wy),
        mean = mean(y[ry], na.rm = TRUE),
        sd = fallback_sd
      )
    )
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
  
  if (is.na(residual_sd) || residual_sd == 0) {
    residual_sd <- sd(y[ry], na.rm = TRUE)
  }
  
  if (is.na(residual_sd) || residual_sd == 0) {
    residual_sd <- 1
  }
  
  imputed_values <- predicted_values + rnorm(
    n = sum(wy),
    mean = 0,
    sd = residual_sd
  )
  
  return(as.numeric(imputed_values))
}


# -----------------------------------------

# I added this because each outcome type needs a different true beta.
# For the count cases, the beta is on the log scale.

true_beta_for_case <- function(case_type) {
  
  if (case_type == "Gaussian") {
    return(1.8)
  }
  
  if (case_type == "Overdispersed") {
    return(0.08)
  }
  
  if (case_type == "Zero-inflated") {
    return(0.08)
  }
}


# -----------------------------------------

sim_data <- function(nschools,
                     n_per_school,
                     true_beta_study_hours,
                     case_type) {
  
  # I made the dataset larger than before but we can still make it even larger
  # When we make it larger it mice() takes much more time to run the itrerations
  
  school_id <- as.integer(rep(1:nschools, each = n_per_school))
  
  age <- rnorm(nschools * n_per_school, mean = 16, sd = 1)
  study_hours <- rnorm(nschools * n_per_school, mean = 5, sd = 2)
  
  
  if (case_type == "Gaussian") {
    
    # Random school effects
    
    school_intercept <- rnorm(nschools, mean = 0, sd = 4)
    school_slope <- rnorm(nschools, mean = 0, sd = 0.4)
    
    u_school <- school_intercept[school_id]
    u_slope <- school_slope[school_id]
    
    # True score
    score_true <- 50 +
      2.5 * age +
      true_beta_study_hours * study_hours +
      u_school +
      u_slope * study_hours +
      rnorm(nschools * n_per_school, mean = 0, sd = 5)
  }
  
  
  if (case_type == "Overdispersed") {
    
    # This case makes a count outcome instead of a normal continuous outcome.
    # I first make the expected count using a log link, then I draw the score from a negative binomial distribution.
    # The negative binomial distribution makes the data overdispersed because the variance is larger than the mean.
    
    school_intercept <- rnorm(nschools, mean = 0, sd = 0.25)
    school_slope <- rnorm(nschools, mean = 0, sd = 0.03)
    
    u_school <- school_intercept[school_id]
    u_slope <- school_slope[school_id]
    
    eta <- 1 +
      0.03 * age +
      true_beta_study_hours * study_hours +
      u_school +
      u_slope * study_hours
    
    mu <- exp(eta)
    
    score_true <- rnbinom(
      nschools * n_per_school,
      mu = mu,
      size = 2
    )
  }
  
  
  if (case_type == "Zero-inflated") {
    
    # This case also makes a count outcome.
    # I first make an overdispersed count outcome like above.
    # Then I randomly force some observations to be zero to create extra zeros in the dataset.
    
    school_intercept <- rnorm(nschools, mean = 0, sd = 0.25)
    school_slope <- rnorm(nschools, mean = 0, sd = 0.03)
    
    u_school <- school_intercept[school_id]
    u_slope <- school_slope[school_id]
    
    eta <- 1 +
      0.03 * age +
      true_beta_study_hours * study_hours +
      u_school +
      u_slope * study_hours
    
    mu <- exp(eta)
    
    count_score <- rnbinom(
      nschools * n_per_school,
      mu = mu,
      size = 2
    )
    
    extra_zero <- rbinom(
      nschools * n_per_school,
      size = 1,
      prob = 0.35
    )
    
    score_true <- ifelse(extra_zero == 1, 0, count_score)
  }
  
  
  data.frame(
    school_id = school_id,
    age = age,
    study_hours = study_hours,
    score = score_true
  )
}


# -----------------------------------------

get_estimates <- function(model, true_beta) {
  
  if (inherits(model, "glmmTMB")) {
    
    coef_table <- summary(model)$coefficients$cond
    
    if (!"study_hours" %in% rownames(coef_table)) {
      return(
        c(
          estimate = NA,
          se = NA,
          covered = NA,
          ciw = NA
        )
      )
    }
    
    est <- coef_table["study_hours", "Estimate"]
    se <- coef_table["study_hours", "Std. Error"]
    
    cov <- est - 1.96 * se < true_beta & true_beta < est + 1.96 * se
    
    return(
      c(
        estimate = est,
        se = se,
        covered = cov,
        ciw = 2 * 1.96 * se
      )
    )
  }
}


# -----------------------------------------

fit_model <- function(data, true_beta, case_type) {
  
  if (case_type == "Gaussian") {
    
    fit <- try(
      glmmTMB::glmmTMB(
        score ~ age + study_hours + (1 + study_hours | school_id),
        data = data,
        REML = FALSE
      ),
      silent = TRUE
    )
  }
  
  
  if (case_type == "Overdispersed") {
    
    # For this case I fit the analysis model as a negative binomial mixed model.
    # This matches the way I generated the overdispersed count outcome.
    
    fit <- try(
      glmmTMB::glmmTMB(
        score ~ age + study_hours + (1 + study_hours | school_id),
        data = data,
        family = glmmTMB::nbinom2
      ),
      silent = TRUE
    )
  }
  
  
  if (case_type == "Zero-inflated") {
    
    # For this case I fit the model with ziformula = ~ 1.
    # This tells glmmTMB to estimate one extra-zero probability for the outcome.
    
    fit <- try(
      glmmTMB::glmmTMB(
        score ~ age + study_hours + (1 + study_hours | school_id),
        ziformula = ~ 1,
        data = data,
        family = glmmTMB::nbinom2
      ),
      silent = TRUE
    )
  }
  
  
  if (inherits(fit, "try-error")) {
    return(
      c(
        estimate = NA,
        se = NA,
        covered = NA,
        ciw = NA
      )
    )
  }
  
  get_estimates(fit, true_beta)
}


# -----------------------------------------

# low priority: argument to use lme4 for fitting lme4 imp case (appendix)

imp <- function(miss_data,
                true_beta,
                method,
                predmat,
                case_type,
                m = 5,
                maxit = 5) {
  
  imp <- mice::mice(
    miss_data,
    method = method,
    predictorMatrix = predmat,
    m = m,
    maxit = maxit,
    printFlag = FALSE
  )
  
  comp <- mice::complete(imp, "all")
  
  out <- sapply(
    comp,
    \(x) fit_model(
      data = x,
      true_beta = true_beta,
      case_type = case_type
    )
  )
  
  ok <- !is.na(out["estimate", ]) & !is.na(out["se", ])
  
  if (sum(ok) < 2) {
    return(
      c(
        estimate = NA,
        se = NA,
        covered = NA,
        ciw = NA
      )
    )
  }
  
  pooled <- mice::pool.scalar(
    Q = out["estimate", ok],
    U = out["se", ok]^2
  )
  
  # hypotheiss test on multilevel data (t dist on reg coeff: check df)
  # The degrees of freedom come from mice::pool.scalar after pooling the estimates
  tquantile <- qt(0.975, pooled$df) * sqrt(pooled$t)
  
  c(
    estimate = pooled$qbar,
    se = sqrt(pooled$t),
    covered = pooled$qbar - tquantile < true_beta & true_beta < pooled$qbar + tquantile,
    ciw = 2 * tquantile
  )
}


# -----------------------------------------

cl <- parallel::makeCluster(16)

parallel::clusterEvalQ(cl, {
  library(glmmTMB)
  library(mice)
  library(miceadds)
  library(lme4)
})

parallel::clusterExport(
  cl,
  varlist = c(
    "mice.impute.2l.glmmTMB",
    "true_beta_for_case",
    "sim_data",
    "get_estimates",
    "fit_model",
    "imp"
  )
)

nsim <- 100

# I added case_type here so each simulation is repeated for the normal, overdispersed, and zero-inflated cases.

case_types <- c(
  "Gaussian",
  "Overdispersed",
  "Zero-inflated"
)

sim_grid <- expand.grid(
  sim = 1:nsim,
  case_type = case_types,
  stringsAsFactors = FALSE
)

sim_out <- pblapply(
  1:nrow(sim_grid),
  \(i) {
    
    sim_id <- sim_grid$sim[i]
    case_type <- sim_grid$case_type[i]
    
    set.seed(2026 + sim_id + match(case_type, case_types) * 100000)
    
    # True parameter I care about: effect of study_hours on score
    true_beta_study_hours <- true_beta_for_case(case_type)
    
    comp <- sim_data(
      100,
      50,
      true_beta_study_hours,
      case_type
    )
    
    comp_est <- fit_model(
      comp,
      true_beta_study_hours,
      case_type
    )
    
    miss <- comp
    
    # I keep score complete because the added count cases would need a separate count imputation model.
    # So here the comparison is about imputing missing predictors in the different outcome settings.
    
    miss[sample(1:nrow(comp), 0.20 * nrow(comp)), "age"] <- NA
    miss[sample(1:nrow(comp), 0.25 * nrow(comp)), "study_hours"] <- NA
    
    cca_est <- fit_model(
      miss,
      true_beta_study_hours,
      case_type
    )
    
    method_glmmTMB <- c(
      school_id = "",
      age = "2l.glmmTMB",
      study_hours = "2l.glmmTMB",
      score = ""
    )
    
    method_lme4 <- c(
      school_id = "",
      age = "2l.lmer",
      study_hours = "2l.lmer",
      score = ""
    )
    
    pred <- mice::make.predictorMatrix(miss)
    pred[, "school_id"] <- -2
    
    # I allow random slopes in the imputation models for the missing predictors.
    pred["age", "study_hours"] <- 2
    pred["study_hours", "age"] <- 2
    
    run_time_glmmtmb <- system.time({
      glmmTMB_est <- imp(
        miss,
        true_beta_study_hours,
        method = method_glmmTMB,
        predmat = pred,
        case_type = case_type
      )
    })
    
    run_time_lmer <- system.time({
      lmer_est <- imp(
        miss,
        true_beta_study_hours,
        method = method_lme4,
        predmat = pred,
        case_type = case_type
      )
    })
    
    
    out <- data.frame(
      case_type = case_type,
      method = c(
        "Complete data",
        "Complete-case analysis",
        "glmmTMB",
        "lmer"
      ),
      rbind(
        t(comp_est),
        t(cca_est),
        t(glmmTMB_est),
        t(lmer_est)
      ),
      row.names = NULL
    )
    
    out$bias <- out$estimate - true_beta_study_hours
    
    out$runtime <- c(
      NA,
      NA,
      run_time_glmmtmb["elapsed"],
      run_time_lmer["elapsed"]
    )
    
    out
  },
  cl = cl
)

parallel::stopCluster(cl)

saveRDS(sim_out, file = "C:/Users/5868777/Downloads/sim_badr_20260625.rds")

sim_out <- do.call(rbind, sim_out)

sim_out$absolute_bias <- abs(sim_out$bias)

simulation_summary <- aggregate(
  cbind(estimate, se, covered, ciw, bias, absolute_bias) ~ case_type + method,
  data = sim_out,
  FUN = function(x) mean(x, na.rm = TRUE)
)

# Runtime is separate because Full data and CCA have NA runtime.
runtime_summary <- aggregate(
  runtime ~ case_type + method,
  data = sim_out,
  FUN = function(x) mean(x, na.rm = TRUE)
)

names(runtime_summary)[names(runtime_summary) == "runtime"] <- "runtime_seconds"

simulation_summary <- merge(
  simulation_summary,
  runtime_summary,
  by = c("case_type", "method"),
  all.x = TRUE
)

simulation_summary


# -----------------------------------------
# Export graphs for report

graph_folder <- "C:/Users/5868777/Downloads/simulation_graphs"

if (!dir.exists(graph_folder)) {
  dir.create(graph_folder)
}

# I keep the method order fixed so the graphs are easier to compare.

simulation_summary$method <- factor(
  simulation_summary$method,
  levels = c(
    "Complete data",
    "Complete-case analysis",
    "glmmTMB",
    "lmer"
  )
)

simulation_summary$case_type <- factor(
  simulation_summary$case_type,
  levels = c(
    "Gaussian",
    "Overdispersed",
    "Zero-inflated"
  )
)


# -----------------------------------------
# Graph 1: Mean absolute bias

# Smaller absolute bias means the estimate is closer to the true value.
# I use facets so all outcome cases are in one clean report graph.

graph_bias <- ggplot(
  simulation_summary,
  aes(x = method, y = absolute_bias)
) +
  geom_col(width = 0.65) +
  facet_wrap(~ case_type, scales = "free_y") +
  labs(
    title = "Mean absolute bias by method and outcome case",
    x = NULL,
    y = "Mean absolute bias"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

graph_bias

ggsave(
  filename = file.path(graph_folder, "graph_1_mean_absolute_bias.png"),
  plot = graph_bias,
  width = 8,
  height = 5,
  dpi = 300
)

# -----------------------------------------
# Graph 2: Coverage rate

# This is the proportion of simulation repetitions where the 95% confidence interval
# included the true value.

graph_coverage <- ggplot(
  simulation_summary,
  aes(x = method, y = covered)
) +
  geom_col(width = 0.65) +
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed"
  ) +
  facet_wrap(~ case_type) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    title = "Coverage rate by method and outcome case",
    x = NULL,
    y = "Coverage rate"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

graph_coverage

ggsave(
  filename = file.path(graph_folder, "graph_2_coverage_rate.png"),
  plot = graph_coverage,
  width = 8,
  height = 5,
  dpi = 300
)

# -----------------------------------------
# Graph 3: Runtime

# I only plot runtime for the two mice imputation methods.
# Full data and CCA do not have mice imputation runtime.

runtime_to_plot <- subset(
  simulation_summary,
  !is.na(runtime_seconds)
)

graph_runtime <- ggplot(
  runtime_to_plot,
  aes(x = method, y = runtime_seconds)
) +
  geom_col(width = 0.65) +
  facet_wrap(~ case_type, scales = "free_y") +
  labs(
    title = "Runtime comparison for imputation methods",
    x = NULL,
    y = "Runtime in seconds"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

graph_runtime

ggsave(
  filename = file.path(graph_folder, "graph_3_runtime_comparison.png"),
  plot = graph_runtime,
  width = 8,
  height = 5,
  dpi = 300
)


# -----------------------------------------
# Export tables too

# I also export the summary table because it will be useful when writing the report.

write.csv(
  simulation_summary,
  file = file.path(graph_folder, "simulation_summary.csv"),
  row.names = FALSE
)

write.csv(
  sim_out,
  file = file.path(graph_folder, "simulation_full_results.csv"),
  row.names = FALSE
)

View(sim_out)
View(simulation_summary)

# End of code

# check REML arg (keep it consistent)
# overdispersed case is included now using a negative binomial outcome
# zero-inflated case is included now by forcing extra zeros and fitting ziformula = ~ 1
# hypotheiss test on multilevel data (t dist on reg coeff: check df)
# final concise function
# if you cant do linear algebra directly in R just do it in C++ and call it in R (esp for nested loops)
