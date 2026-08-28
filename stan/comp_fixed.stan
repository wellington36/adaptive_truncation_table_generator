functions {
    real logFunction(real loglamb, real nu, int n) {
        return (n - 1) * loglamb - nu * lgamma(n);
    }

    array[] real sequential(real loglamb, real nu, int MAX_ITERS) {
        vector[MAX_ITERS] storeVal;
        real logZ;
        int i = 1;
  
          storeVal[1] = - nu * lgamma(1);
  
        while (i < MAX_ITERS) {
            i+=1;
            storeVal[i] = logFunction(loglamb, nu, i);
        }
  
        logZ = log_sum_exp(sort_asc(storeVal[:i]));
        return {logZ, 1. * i};
    }
}

data {
  int<lower=1> N;                // Number of unique counts
  array[N] int<lower=0> y;       // Counts (0, 1, 2, ...)
  array[N] int<lower=0> freq;    // Frequencies for each count
  int<lower=1> FIXED;
}

parameters {
  real mu;              // Mean parameter (mu)
  real<lower=0> nu;              // Dispersion parameter (nu)
}

transformed parameters {
  real logZ;                        // Normalization constant
  real loglamb = nu * log(mu);
  array[2] real infiniteSequentialApproach = sequential(loglamb, nu, FIXED);
  
  logZ = infiniteSequentialApproach[1];
}

model {
  vector[N] log_p;          // Log probabilities for each count
  // Priors (adjust these based on your knowledge)
  mu ~ gamma(1, 1);            // Prior for mu
  nu ~ gamma(0.0625, 0.25);        // Prior for nu
  //nu ~ normal(0, 1);               // Prior for nu
  
  // Compute log probabilities
  for (i in 1:N) {
    log_p[i] = y[i] * loglamb - nu * lgamma(y[i] + 1) - logZ;
    target += freq[i] * log_p[i];
  }
}
