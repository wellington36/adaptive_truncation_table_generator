from mpmath import mpf, log, expm1
from utils.utils import logsumexp, logdiffexp

def integration_bound_mp(f, g, theta, M, eps, initial_k):
    k = initial_k
    leps = log(eps)
    log_terms = [log(mpf(0))] * (M + initial_k)

    if (g(theta, M) > log(mpf(2)) + leps):
        raise ValueError("It is not possible to reach the stopping criterion with the given M.")

    log_terms[k] = f(theta, k)
    k+=1

    while (g(theta, k) > log(mpf(2)) + leps):
        log_terms[k] = f(theta, k)
        k+=1
    
    log_sum = logsumexp(log_terms)

    return (k-initial_k, logsumexp([log_sum, logsumexp([g(theta, k), g(theta, k-1)]) - log(mpf(2))]))

