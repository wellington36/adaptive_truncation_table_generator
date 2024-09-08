functions {
  real f(real a, int n) {
    return - 2 * log(n+1) -  (n+1)*log(a);
  }
  int normalising_f_btt(real a, real eps) {
    real L = 1/a;
    real M = (L+1)/2;
    real leps = log(eps);
    int MAX_TERMS = 1000000;
    int k = 1;
    vector[MAX_TERMS] log_terms;
    
    log_terms[1] = f(a, k);
    while (f(a, k) >= leps + log(1-M) - log(M) && k < MAX_TERMS) {
      k += 1;
      log_terms[k] = f(a, k);
    }
    
    return k;
  }
  int normalising_f_bp(real a, real eps) {
    real L = 1/a;
    real M = (L+1)/2;
    real leps = log(eps);
    int MAX_TERMS = 1000000;
    int k = 2;
    vector[MAX_TERMS] log_terms;
    
    log_terms[1] = f(a, 1);
    log_terms[2] = f(a, 2);
    while (log(exp(log_terms[k]) * abs((1-L)^-1 - (1 - exp(log_terms[k])/exp(log_terms[k-1]))^-1)) >= log(2) + leps && k < MAX_TERMS) {
      k += 1;
      log_terms[k] = f(a, k);
    }
    
    return k;
  }
}