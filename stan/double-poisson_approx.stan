functions {
  real log_kernel(real mu, real phi, int n) {
    if (n == 0)
      return 0;
    return -n
            + n * log(n)
            - lgamma(n + 1)
            + phi * n
            + phi * n * log(mu)
            - phi * n * log(n);
  }
}

data {
  int<lower=1> N;
  array[N] int<lower=0> y;
  array[N] int<lower=0> freq;
}

parameters {
  real<lower=0> mu;
  real<lower=0> phi;
}

model {
  vector[N] log_p;
  mu  ~ normal(0, 5);     // constraint handles positivity; fine as a weakly-informative prior
  phi ~ uniform(0, 10);      // replaces the hard-truncated uniform(0,10)

  for (i in 1:N) {
    log_p[i] = 0.5 * log(phi) - phi * mu + log_kernel(mu, phi, y[i]);
    target += freq[i] * log_p[i];
  }
}
