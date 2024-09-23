from mpmath import mp, mpf, exp, log, expm1, fabs
from utils.utils import logsumexp
from math import lgamma
import mpmath

def brute_mp(f, theta, M, initial_k):
    k = initial_k
    log_terms = [log(mpf(0))] * (M+initial_k+1)

    log_terms[k] = f(theta, k)
    log_terms[k+1] = f(theta, k+1)
    k+=1

    while (k < M+initial_k):
        k+=1
        log_terms[k] = f(theta, k)
    
    log_sum = logsumexp(log_terms[initial_k:(k+1)])

    return (k-initial_k, log_sum)