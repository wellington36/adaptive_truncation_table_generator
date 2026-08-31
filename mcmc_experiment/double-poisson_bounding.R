library(rstan)
library(readr)
library(dplyr)
library(knitr)
library(kableExtra)

# MCMC settings
iterations <- 5000
max_iters <- 10^4
eps <- 2^(-52)

# Set rstan options for better performance
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Read your CSV data
data_df <- read_csv("data/Shmuelli_2005.csv")

# Inspect the data
print(data_df)

# Prepare data for Stan
data_df <- data_df %>%
  arrange(count)

counts <- data_df$count
frequencies <- data_df$frequency

# Prepare the list of data for Stan
stan_data <- list(
  N = length(counts),
  y = counts,
  freq = frequencies,
  eps = eps,
  MAX_ITERS = max_iters
)

# Compile the Stan model
stan_model <- stan_model(file = "stan/double-poisson_bounding.stan")

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

# Print a summary of the results
parameters <- c("mu", "phi", "n")
print(fit, pars = parameters)

summary_fit <- summary(fit, pars = parameters)
posterior_stats <- as.data.frame(summary_fit$summary)

# Get elapsed time for each chain
chain_times <- get_elapsed_time(fit)

# Calculate the average time in minutes across all chains
avg_time_min <- mean(rowSums(chain_times)) / 60

# Calculate ESS/minute by dividing n_eff by the average time in minutes
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
