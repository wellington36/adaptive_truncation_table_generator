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

    array[] real infiniteBoundingPairs(real mu, real phi, real epsilon, int MAX_ITERS) {
        vector[MAX_ITERS] storeVal;
        real leps = log(epsilon);
        real logZ;
        int i = 1;
  
        storeVal[1] = log_kernel(mu, phi, 0);
        i+=1;
        storeVal[i] = log_kernel(mu, phi, i-1);
  
        while ((storeVal[i] >= storeVal[i-1] || (storeVal[i] - log(-expm1(storeVal[i] - storeVal[i-1])) >= leps)) && (i < MAX_ITERS)) {
            i+=1;
            storeVal[i] = log_kernel(mu, phi, i-1);
        }
  
        logZ = log_sum_exp(sort_asc(storeVal[:i]));
        return {logZ, 1. * i};
    }
}

data {
  int<lower=1> N;                // Number of unique counts
  array[N] int<lower=0> y;       // Counts (0, 1, 2, ...)
  array[N] int<lower=0> freq;    // Frequencies for each count
  real eps;
  int MAX_ITERS;
}

parameters {
  real<lower=0> mu;              // Mean parameter (mu)
  real<lower=0> phi;              // Dispersion parameter (phi)
}

transformed parameters {
  real logZ;                        // Normalization constant
  array[2] real infiniteBPApproach = infiniteBoundingPairs(mu, phi, eps, MAX_ITERS);
  
  logZ = infiniteBPApproach[1];
}

model {
  vector[N] log_p;
  mu ~ normal(0, 5);
  phi ~ uniform(0, 10);

  for (j in 1:N) {
    log_p[j] = log_kernel(mu, phi, y[j]) - logZ;
    target += freq[j] * log_p[j];
  }
}

generated quantities{
  real n = infiniteBPApproach[2];
}
