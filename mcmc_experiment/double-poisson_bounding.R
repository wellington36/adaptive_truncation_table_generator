library(rstan)
library(readr)
library(dplyr)
library(ggplot2)

iterations = 5000
eps = 2**(-52)

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
  MAX_ITERS = 10**4
)

# Compile the Stan model
stan_model <- stan_model(file = "stan/double-poisson_bounding.stan")

# Fit the model using MCMC
fit <- sampling(
  object = stan_model,
  data = stan_data,
  refresh = floor(iterations/5),
  iter = 2*iterations,               # Number of iterations
  warmup = iterations,    # Number of warmup (burn-in) iterations
  chains = 4,                      # Number of chains
  cores = 8,
  control = list(adapt_delta = 0.90, max_treedepth = 12)  # Control parameters
)

# Print a summary of the results
print(fit, pars = c("mu", "phi", "n"))

summary_fit <- summary(fit, pars = c("mu", "phi", "n"))

# Convert the summary output to a data frame
posterior_stats <- as.data.frame(summary_fit$summary)

# Get elapsed time for each chain
chain_times <- get_elapsed_time(fit)

# Calculate the average time in minutes across all chains
avg_time_min <- mean(rowSums(chain_times)) / 60

# Calculate ESS/minute by dividing n_eff by the average time in minutes
ess_per_minute <- posterior_stats$n_eff / avg_time_min


# Create a summary table for mu and phi
summary_table <- data.frame(
  Parameter = c("mu", "phi", "n"),
  Mean = posterior_stats$mean,
  Median = posterior_stats$`50%`,
  `95% BCI` = paste0("[", round(posterior_stats$`2.5%`, 3), ", ", round(posterior_stats$`97.5%`, 3), "]"),
  `Posterior SD` = posterior_stats$sd,
  MCSE = posterior_stats$se_mean,
  `ESS/minute` = ess_per_minute
)

# Display the summary table
print(summary_table)

# Optional: Format the table for display
library(knitr)
library(kableExtra)

summary_table %>%
  kable("html", col.names = c("Parameter", "Mean", "Median", "95% BCI", "Posterior SD", "MCSE", "ESS/minute")) %>%
  kable_styling(full_width = F, bootstrap_options = c("striped", "hover"))
