from mpmath import mp, mpf, exp, log, expm1
from utils.utils import logsumexp
from math import lgamma
import mpmath

def sum_to_threshold_mp(f, theta, M, L, eps, initial_k):
    k = initial_k
    M_bound = (L+1)/2
    leps = log(eps)
    log_terms = [log(mpf(0))] * (M+initial_k)

    log_terms[k] = f(theta, k)
    log_terms[k+1] = f(theta, k+1)
    k+=1

    while (log_terms[k] > leps + log(1-M_bound) - log(M_bound) and k < M-1+initial_k):
        k+=1
        log_terms[k] = f(theta, k)
    
    log_sum = logsumexp(log_terms[initial_k:(k+1)])

    return (k, log_sum)