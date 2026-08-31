functions {
  // Log of the unnormalized COM-Poisson mass at count y.
  real log_mass(real log_lambda, real nu, int y) {
    return y * log_lambda - nu * lgamma(y + 1);
  }

  // Calculate log(Z) by summing COM-Poisson mass terms in log space.
  // The second returned value is the number of terms used.
  array[] real infiniteBoundingPairs(
    real log_lambda,
    real nu,
    real epsilon,
    int MAX_ITERS
  ) {
    vector[MAX_ITERS] log_mass_values;
    real log_epsilon = log(epsilon);
    real log_Z;
    int i = 1;

    // The first two terms correspond to counts 0 and 1.
    log_mass_values[i] = log_mass(log_lambda, nu, 0);
    i += 1;
    log_mass_values[i] = log_mass(log_lambda, nu, 1);

    // Continue until the next term is decreasing and negligible relative
    // to the previous term, or until MAX_ITERS is reached.
    while (
      (
        log_mass_values[i] >= log_mass_values[i - 1] ||
        log_mass_values[i] - log(-expm1(log_mass_values[i] - log_mass_values[i - 1])) >= log_epsilon
      ) && i < MAX_ITERS
    ) {
      i += 1;
      log_mass_values[i] = log_mass(log_lambda, nu, i - 1);
    }

    log_Z = log_sum_exp(sort_asc(log_mass_values[:i]));
    return {log_Z, 1.0 * i};
  }
}

data {
  int<lower=1> N;             // Number of unique observed counts.
  array[N] int<lower=0> y;    // Observed count values.
  array[N] int<lower=0> freq; // Frequency associated with each count.
  real<lower=0> eps;          // Truncation tolerance.
  int<lower=2> MAX_ITERS;     // Maximum number of terms in Z.
}

parameters {
  real<lower=0> mu;           // COM-Poisson location parameter.
  real<lower=0> nu;           // COM-Poisson dispersion parameter.
  real<lower=0, upper=1> zi;  // Probability of an additional structural zero.
}

transformed parameters {
  real log_lambda;
  real log_Z;
  real n;
  array[2] real bounding_pair;

  log_lambda = nu * log(mu);
  bounding_pair = infiniteBoundingPairs(log_lambda, nu, eps, MAX_ITERS);
  log_Z = bounding_pair[1];
  n = bounding_pair[2];
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

generated quantities {
  // Number of terms used to approximate the COM-Poisson normalizing constant.
  real truncation_terms = n;
}
