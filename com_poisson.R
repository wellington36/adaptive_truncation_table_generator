library(rstan)
library(COMPoissonReg)

expose_stan_functions("stan/comp_normalization_bounding.stan")
expose_stan_functions("stan/com_poisson_brute.stan")

mu <- c(10, 100, 1000, 10000)
nu <- c(0.1, 0.01, 0.001, 0.0001)
M <- c(10^3, 10^4, 10^5, 3*10^5) # only for the brute (10^6 for bounding)
lambda <- mu^(nu)

brute_lessbits_error <- numeric(length(lambda))
for (i in seq_along(lambda)) {
  brute_lessbits_error[i] <- log_Z_brute(log(lambda[i]), nu[i], M[i], 2^-52*10^6)
}

bounding_lessbits_error <- numeric(length(lambda))
for (i in seq_along(lambda)) {
  bounding_lessbits_error[i] <- log_Z_bounding(log(lambda[i]), nu[i], 10^6, 2^-52*10^6)
}

brute_64bits_error <- numeric(length(lambda))
for (i in seq_along(lambda)) {
  brute_64bits_error[i] <- log_Z_brute(log(lambda[i]), nu[i], M[i], 2^-52)
}

bounding_64bits_error <- numeric(length(lambda))
for (i in seq_along(lambda)) {
  bounding_64bits_error[i] <- log_Z_bounding(log(lambda[i]), nu[i], 10^6, 2^-52)
}


table <- data.frame(
  Column1 = brute_lessbits_error,
  Column2 = bounding_lessbits_error,
  Column3 = brute_64bits_error,
  Column4 = bounding_64bits_error
)

# Set the column names
colnames(table) <- c("[2.2x10^-10 | Brute]", "[2.2x10^-10 | Bounding Pairs]", "[2.2x10^-16 | Brute]", "[2.2x10^-16 | Bounding Pairs]")

# Set the row names
rownames(table) <- c("mu=10 | nu=0.1", "mu=100 | nu=0.01", "mu=1000 | nu=0.001", "mu=10000 | nu=0.0001")

# Print the table
print(table)