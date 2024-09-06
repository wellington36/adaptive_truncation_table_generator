# code by: https://github.com/paul-buerkner/brms/issues/607
library(rstan)
library(COMPoissonReg)

#expose_stan_functions("com_poisson_brms.stan")
#expose_stan_functions("com_poisson_max.stan")
#expose_stan_functions("com_poisson_brute.stan")
expose_stan_functions("stan/com_poisson_bounding.stan")
#expose_stan_functions("com_poisson_minka.stan")
#expose_stan_functions("com_poisson_true.stan")
#expose_stan_functions("sandbox.stan")
expose_stan_functions("stan/comp_normalization_bounding.stan")

y <- 0
mu <- 10
nu <- 0.1
lambda <- mu^(nu)

#br <- numeric(length(y))

#M <- 100
#print(M * log(mu) - nu*lgamma(M + 1))
#print((mu**M)/((factorial(M)**nu)))
#print(f(log(lambda^(1/nu)), nu))


#for (i in seq_along(y)) {
#  br[i] <- exp(com_poisson_lpmf(y[i], mu, nu))
#}

#br_max <- numeric(length(y))
#for (i in seq_along(y)) {
# br_max[i] <- exp(com_poisson_lpmf_max(y[i], lambda, nu))
#}

#br_bounding <- numeric(length(y))
#for (i in seq_along(y)) {
# br_bounding[i] <- exp(com_poisson_lpmf_bounding(y[i], lambda, nu))
#}

bounding <- numeric(length(y))
for (i in seq_along(y)) {
  bounding[i] <- exp(log_Z_com_poisson_bounding(log(lambda), nu))
}

b <- log_Z_bounding(log(lambda), nu, 10000, 2^-52)
print(b)

#br_brute <- numeric(length(y))
#for (i in seq_along(y)) {
#  br_brute[i] <- exp(com_poisson_lpmf_brute(y[i], mu, nu))
#}

#br_true <- numeric(length(y))
#for (i in seq_along(y)) {
# br_true[i] <- exp(com_poisson_lpmf_true(y[i], lambda, nu))
#}

#br_minka <- numeric(length(y))
#for (i in seq_along(y)) {
#  br_minka[i] <- exp(com_poisson_lpmf_minka(y[i], lambda, nu))
#}

results <- data.frame(
  y = y,
  #brms = br,
  #brms_max = br_max,
  brms_bounding = bounding
  #brms_true = br_true,
  #diff = br_max - br_bounding,
  #brms_minka = br_minka,
  #COMPoissonReg = dcmp(y, lambda, nu),
  #error_minka = abs(br_brute - br_minka),
  #error_COMPoissonReg = abs(br_brute - dcmp(y, lambda, nu))
)

print(results)
