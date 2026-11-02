source("experiments.R")

#canada_results <- run_published_data(
#  "blumberg_canada_chain_sizes.csv",
#  size_column = "chain_size"
#)

us_results <- run_published_data(
  "blumberg_us_chain_sizes.csv",
  size_column = "chain_size",
  fixed_J = c(10L, 50L, 100L, 500L),
  eps_ll = 1e-10
)


head(us_results)
plot_case_study(us_results)