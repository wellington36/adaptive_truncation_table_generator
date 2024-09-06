library(rstan)
library(COMPoissonReg)

#expose_stan_functions("stan/comp_normalization_bounding.stan")

mu <- c(10, 100, 1000, 10000)
nu <- c(0.1, 0.01, 0.001, 0.0001)
lambda <- mu^(nu)

print(log(lambda))

bounding_64bits_error <- numeric(length(lambda))
for (i in seq_along(lambda)) {
  bounding_64bits_error[i] <- log_Z_bounding(log(lambda[i]), nu[i], 300000, 2^-52)
}

bounding_lessbits_error <- numeric(length(lambda))
for (i in seq_along(lambda)) {
  bounding_lessbits_error[i] <- log_Z_bounding(log(lambda[i]), nu[i], 300000, 2^-52*10^6)
}


table <- data.frame(Column1 = bounding_lessbits_error, Column2 = bounding_64bits_error)

# Set the column names
colnames(table) <- c("[2.2x10^-10 | Bounding Pairs]", "[2.2x10^-16 | Bounding Pairs]")

# Set the row names
rownames(table) <- c("mu=10 | nu=0.1", "mu=100 | nu=0.01", "mu=1000 | nu=0.001", "mu=10000 | nu=0.0001")

# Print the table
print(table)