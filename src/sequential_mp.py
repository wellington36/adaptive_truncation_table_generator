from mpmath import mp, mpf, exp, log, expm1, fabs
from utils.utils import logsumexp, logdiffexp
import mpmath

def sequential_mp(f, theta, M, eps, initial_k):
    k = initial_k
    leps = log(eps)
    log_terms_brute = [log(mpf(0))] * (M+initial_k+1)

    log_terms_brute[k] = f(theta, k)
    log_terms_brute[k+1] = f(theta, k+1)
    k+=1

    while (k < M+initial_k):
        k+=1
        log_terms_brute[k] = f(theta, k)
    
    log_sum_brute = logsumexp(log_terms_brute[initial_k:(k+1)])

    k = initial_k
    log_terms = [log(mpf(0))] * (M+initial_k+1)

    log_terms[k] = f(theta, k)
    log_terms[k+1] = f(theta, k+1)
    k+=1
    log_sum = logsumexp(log_terms[initial_k:(k+1)])

    while (logdiffexp(log_sum_brute, log_sum) >= leps):
        k+=1
        log_terms[k] = f(theta, k)
        log_sum = logsumexp([log_sum, log_terms[k]])
    
    
    return (k-initial_k, log_sum)