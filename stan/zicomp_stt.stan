functions {
    // Log of the unnormalized COM-Poisson mass at count y.
    real log_kernel(real log_lambda, real nu, int y) {
        return y * log_lambda - nu * lgamma(y + 1);
    }
    
    array[] real infiniteSumToThreshold(real log_lambda, real nu, real L,
                                        real epsilon, int M, int initial_k) {
        vector[M] storeVal;          // at most M terms will ever be stored
        real leps = log(epsilon);
        real M_bound = (L + 1) / 2;
        real log_M_bound = log(M_bound);
        real threshold = leps + log1m(M_bound) - log(M_bound);
        real logZ;
        int idx = 1;
        int k = initial_k;
        int ratio_ok = 1;

        if (log_kernel(log_lambda, nu, M) > threshold) {
            reject("It is not possible to reach the stopping criterion with the given M.");
        }

        storeVal[1] = log_kernel(log_lambda, nu, initial_k);

        // The stopping criterion (below threshold) is only valid when the ratio
        // between consecutive terms doesn't exceed M_bound. Track both, and only
        // stop the loop once BOTH conditions hold.
        while ((storeVal[idx] > threshold || ratio_ok == 0) && k < M - 1 + initial_k) {
            k += 1;
            idx += 1;
            storeVal[idx] = log_kernel(log_lambda, nu, k);
            ratio_ok = (storeVal[idx] - storeVal[idx - 1]) <= log_M_bound ? 1 : 0;
        }

        logZ = log_sum_exp(sort_asc(storeVal[:idx]));
        return {logZ, 1. * (k - initial_k)};
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
  bounding_pair = infiniteSumToThreshold(log_lambda, nu, 0, eps, MAX_ITERS, 0);
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
