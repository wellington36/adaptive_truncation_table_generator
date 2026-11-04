source("experiments.R")

canada_results <- run_published_data(
  "blumberg_canada_chain_sizes.csv",
  size_column = "chain_size",
  fixed_J = c(20L, 50L, 200L, 1000L),
  eps_ll = 1e-10
)

us_results <- run_published_data(
  "blumberg_us_chain_sizes.csv",
  size_column = "chain_size",
  fixed_J = c(20L, 50L, 200L, 1000L),
  eps_ll = 1e-10
)


head(us_results)
plot_case_study(us_results)