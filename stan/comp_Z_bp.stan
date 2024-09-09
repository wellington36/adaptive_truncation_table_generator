functions{
  int log_Z_iterations_bp(real log_lambda, real nu, int M, real eps) {
    real log_Z;
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
    while (log(exp(log_Z_terms[k]) * abs(1 - pow(1 - exp(log_Z_terms[k])/exp(log_Z_terms[k-1]), -1))) >= log(2) + leps && k < M) {
      k += 1;
      log_Z_terms[k] = (k - 1) * log_lambda - nu * lgamma(k);
    }

    return k-1;
  }
  
  real log_Z_value_bp(real log_lambda, real nu, int M, real eps) {
    real log_Z;
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
    while (log(exp(log_Z_terms[k]) * abs(1 - pow(1 - exp(log_Z_terms[k] - log_Z_terms[k-1]), -1))) >= log(2) + leps && k < M) {
      k += 1;
      log_Z_terms[k] = (k - 1) * log_lambda - nu * lgamma(k);
    }

    log_Z = log_sum_exp(log_Z_terms[1:k]);
  
    return log_sum_exp([log_Z, log(exp(log_Z_terms[k])/2 * (1 + pow(1 - exp(log_Z_terms[k] - log_Z_terms[k-1]), -1)))]);
  }
}