library(rstan)
library(COMPoissonReg)
source("aux/aux.R")

expose_stan_functions("stan/comp_Z_bp.stan")
expose_stan_functions("stan/comp_Z_brute.stan")
expose_stan_functions("stan/comp_Z_brms_fixed.stan")

mu <- c(10, 100, 1000, 10000)
nu <- c(0.1, 0.01, 0.001, 0.0001)
lambda <- mu^(nu)

# M only for the brute (10^6 for bounding)
M <- c(10^4, 10^4, 10^5, 3*10^5) # test
#M <- c(10^4, 10^5, 10^6, 10^6)  # medium
#M <- c(10^6, 10^6, 10^6, 10^6)  # final version


brute_lessbits_error <- numeric(length(lambda))
for (i in seq_along(lambda)) {
  brute_lessbits_error[i] <- log_Z_iterations_brute(log(lambda[i]), nu[i], M[i], 2^-52*10^6)
}

bounding_lessbits_error <- numeric(length(lambda))
for (i in seq_along(lambda)) {
  bounding_lessbits_error[i] <- log_Z_iterations_bp(log(lambda[i]), nu[i], 10^6, 2^-52*10^6)
}

brute_64bits_error <- numeric(length(lambda))
for (i in seq_along(lambda)) {
  brute_64bits_error[i] <- log_Z_iterations_brute(log(lambda[i]), nu[i], M[i], 2^-52)
}

bounding_64bits_error <- numeric(length(lambda))
for (i in seq_along(lambda)) {
  bounding_64bits_error[i] <- log_Z_iterations_bp(log(lambda[i]), nu[i], 10^6, 2^-52)
}

library_COMPoissonReg <- numeric(length(lambda))
for (i in seq_along(lambda)) {
  library_COMPoissonReg[i] <- abs(exp(log_diff_exp(log_Z_brute(log(lambda[i]), nu[i], M[i]), -1*dcmp(0, lambda[i], nu[i], log = TRUE))))
  print(exp(-1*dcmp(0, lambda[i], nu[i], log = TRUE)))
  print(exp(log_Z_brute(log(lambda[i]), nu[i], M[i])))
}

library_brms <- numeric(length(lambda))
for (i in seq_along(lambda)) {
  library_brms[i] <- abs(exp(log_diff_exp(log_Z_brute(log(lambda[i]), nu[i], M[i]), log_Z_com_poisson(log(mu[i]), nu[i]))))
}

table <- data.frame(
  Column1 = brute_lessbits_error,
  Column2 = bounding_lessbits_error,
  Column3 = brute_64bits_error,
  Column4 = bounding_64bits_error,
  Column5 = library_brms,
  Column6 = library_COMPoissonReg
)

# Set the column names
colnames(table) <- c("[2.2x10^-10|Brute]", "[2.2x10^-10|BP]", "[2.2x10^-16|Brute]", "[2.2x10^-16|BP]", "brms", "COMPoissonReg")

# Set the row names
rownames(table) <- c("mu=10, nu=0.1", "mu=100, nu=0.01", "mu=1000, nu=0.001", "mu=10000, nu=0.0001")

# Print the table
print(table)