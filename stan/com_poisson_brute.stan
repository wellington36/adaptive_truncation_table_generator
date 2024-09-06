functions {
  int log_Z_brute(real log_lambda, real nu, int M, real eps) {
    real log_Z;
    real log_Z_true;
    int k = 2;
    real leps = log(eps);
    vector[M] log_Z_terms;

    if (nu == 1) {
      reject("nu = 1");
    }
    if (nu <= 0) {
      reject("nu must be positive");
    }
    if (nu == positive_infinity()) {
      reject("nu must be finite");
    }

    log_Z_terms[1] = 0;
    log_Z_terms[2] = log_lambda;
    while (k < M) {
      k += 1;
      log_Z_terms[k] = (k - 1) * log_lambda - nu * lgamma(k);
    }
    log_Z = log_sum_exp(log_Z_terms[1:k]);
    
    vector[M] log_Z_true_terms;
    k = 2;
    log_Z_true_terms[1] = 0;
    log_Z_true_terms[2] = log_lambda;
    log_Z_true = log_sum_exp(log_Z_true_terms[1:k]);
    

    while (log_diff_exp(log_Z, log_Z_true) >= leps) {
      k += 1;
      log_Z_true_terms[k] = (k - 1) * log_lambda - nu * lgamma(k);
      log_Z_true = log_sum_exp(log_Z_true_terms[1:k]);
    }

    return k-1;
  }
}