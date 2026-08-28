library(rstan)
library(COMPoissonReg)
library(dplyr)
# Optional packages used for formatted output.
library(knitr)
library(kableExtra)

# MCMC settings.
iterations <- 5000
max_iters <- 10^4
eps <- 2^(-52)

# Configure rstan for better performance.
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---- Load and prepare the couple$UPB data ----
data(couple)
y_obs <- couple$UPB

# Convert raw counts into grouped count/frequency data, as required
# by the zicompoisson_bounding.stan model.
count_table <- table(y_obs)
data_df <- data.frame(
  count = as.integer(names(count_table)),
  frequency = as.integer(count_table)
) %>%
  arrange(count)

print(data_df)

# Validate the required columns and their values.
stopifnot(all(c("count", "frequency") %in% names(data_df)))
stopifnot(all(data_df$count >= 0))
stopifnot(all(data_df$frequency >= 0))

# Prepare the data for Stan.
counts <- as.integer(data_df$count)
frequencies <- as.integer(data_df$frequency)
stan_data <- list(
  N = length(counts),
  y = counts,
  freq = frequencies,
  eps = eps,
  MAX_ITERS = max_iters
)

# Compile the ZI-COM-Poisson Stan model.
stan_model <- stan_model(file = "stan/zicomp_bounding.stan")

# Fit the model using MCMC.
fit <- sampling(
  object = stan_model,
  data = stan_data,
  refresh = floor(iterations / 5),
  iter = 2 * iterations,
  warmup = iterations,
  chains = 4,
  cores = min(4, parallel::detectCores()),
  control = list(
    adapt_delta = 0.90,
    max_treedepth = 12
  )
)

# Print the posterior summary for all model parameters and generated quantities.
parameters <- c("mu", "nu", "zi", "n")
print(fit, pars = parameters)
summary_fit <- summary(fit, pars = parameters)
posterior_stats <- as.data.frame(summary_fit$summary)

# Calculate the average elapsed time across chains, in minutes.
chain_times <- get_elapsed_time(fit)
avg_time_min <- mean(rowSums(chain_times)) / 60

# Calculate effective sample size per minute.
ess_per_minute <- posterior_stats$n_eff / avg_time_min

# Create a compact posterior summary table.
summary_table <- data.frame(
  Parameter = parameters,
  Mean = posterior_stats$mean,
  Median = posterior_stats$`50%`,
  `95% BCI` = paste0(
    "[",
    round(posterior_stats$`2.5%`, 3),
    ", ",
    round(posterior_stats$`97.5%`, 3),
    "]"
  ),
  `Posterior SD` = posterior_stats$sd,
  MCSE = posterior_stats$se_mean,
  `ESS/minute` = ess_per_minute,
  check.names = FALSE
)
print(summary_table)

# Optional HTML rendering of the summary table.
summary_table %>%
  kable(
    format = "html",
    col.names = c(
      "Parameter", "Mean", "Median", "95% BCI",
      "Posterior SD", "MCSE", "ESS/minute"
    )
  ) %>%
  kable_styling(
    full_width = FALSE,
    bootstrap_options = c("striped", "hover")
  )
