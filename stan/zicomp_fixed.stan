functions {
    // Log of the unnormalized COM-Poisson mass at count y.
    real log_kernel(real log_lambda, real nu, int y) {
        return y * log_lambda - nu * lgamma(y + 1);
    }

    real sequential(real log_lambda, real nu, int max_iters) {
        vector[max_iters + 1] log_terms;

        log_terms[1] = 0;  // n = 0

        for (n in 1:max_iters)
            log_terms[n + 1] = log_kernel(log_lambda, nu, n);

        return log_sum_exp(log_terms);
    }
}

data {
  int<lower=1> N;             // Number of unique observed counts.
  array[N] int<lower=0> y;    // Observed count values.
  array[N] int<lower=0> freq; // Frequency associated with each count.
  int<lower=2> FIXED;     // Maximum number of terms in Z.
}

parameters {
  real<lower=0> mu;           // COM-Poisson location parameter.
  real<lower=0> nu;           // COM-Poisson dispersion parameter.
  real<lower=0, upper=1> zi;  // Probability of an additional structural zero.
}

transformed parameters {
  real log_lambda;
  real log_Z;

  log_lambda = nu * log(mu);
  log_Z = sequential(log_lambda, nu, FIXED);
}

model {
  vector[N] log_p;

  // Priors. Adjust these to reflect substantive prior knowledge.
  mu ~ gamma(1, 1);
  nu ~ gamma(0.0625, 0.25);
  zi ~ beta(1, 1);

  for (j in 1:N) {
    real log_com_poisson_p;

    // Log probability under the ordinary COM-Poisson component.
    log_com_poisson_p =
      y[j] * log_lambda - nu * lgamma(y[j] + 1) - log_Z;

    // A zero may arise either from the structural-zero component or from
    // the ordinary COM-Poisson component. Positive counts can only arise
    // from the ordinary component.
    if (y[j] == 0) {
      log_p[j] = log_mix(zi, 0, log_com_poisson_p);
    } else {
      log_p[j] = log1m(zi) + log_com_poisson_p;
    }

    // Use the grouped frequency as a likelihood multiplier.
    target += freq[j] * log_p[j];
  }
}
