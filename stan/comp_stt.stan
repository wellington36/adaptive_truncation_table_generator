functions {
    // term(k) = k*log(lambda) - nu*lgamma(k+1), i.e. log of lambda^k / (k!)^nu
    real logTerm(real loglamb, real nu, int k) {
        return k * loglamb - nu * lgamma(k + 1);
    }

    array[] real infiniteSumToThreshold(real loglamb, real nu, real L,
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

        if (logTerm(loglamb, nu, M) > threshold) {
            reject("It is not possible to reach the stopping criterion with the given M.");
        }

        storeVal[1] = logTerm(loglamb, nu, initial_k);

        // The stopping criterion (below threshold) is only valid when the ratio
        // between consecutive terms doesn't exceed M_bound. Track both, and only
        // stop the loop once BOTH conditions hold.
        while ((storeVal[idx] > threshold || ratio_ok == 0) && k < M - 1 + initial_k) {
            k += 1;
            idx += 1;
            storeVal[idx] = logTerm(loglamb, nu, k);
            ratio_ok = (storeVal[idx] - storeVal[idx - 1]) <= log_M_bound ? 1 : 0;
        }

        logZ = log_sum_exp(sort_asc(storeVal[:idx]));
        return {logZ, 1. * (k - initial_k)};
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
  real mu;              // Mean parameter (mu)
  real<lower=0> nu;              // Dispersion parameter (nu)
}

transformed parameters {
  real logZ;                        // Normalization constant
  real loglamb = nu * log(mu);
  array[2] real infiniteSTTApproach = infiniteSumToThreshold(loglamb, nu, 0, eps, MAX_ITERS, 0);
  
  logZ = infiniteSTTApproach[1];
}

model {
  vector[N] log_p;          // Log probabilities for each count
  // Priors (adjust these based on your knowledge)
  mu ~ gamma(1, 1);            // Prior for mu
  nu ~ gamma(0.0625, 0.25);        // Prior for nu
  //nu ~ normal(0, 1);               // Prior for nu
  
  // Compute log probabilities
  for (j in 1:N) {
    log_p[j] = y[j] * loglamb - nu * lgamma(y[j] + 1) - logZ;
    target += freq[j] * log_p[j];
  }
}

generated quantities{
  real n = infiniteSTTApproach[2];
}
