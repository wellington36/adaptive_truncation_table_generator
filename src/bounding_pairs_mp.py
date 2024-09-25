from mpmath import mp, mpf, exp, log, expm1
from utils.utils import logsumexp
from math import lgamma
import mpmath

def bounding_pairs_mp(f, theta, M, L, eps, initial_k):
    k = initial_k
    leps = log(eps)
    log_terms = [log(mpf(0))] * (M+initial_k)

    if (f(theta, M) + log(abs(L + (mpf(1)/(mpf(1) - L*(L*L)) + mpf(1)/(expm1(f(theta, M) - f(theta, M-1)))))) >= log(mpf(2)) + leps):
        raise ValueError("It is not possible to reach the stopping criterion with the given M.")

    log_terms[k] = f(theta, k)
    log_terms[k+1] = f(theta, k+1)
    k+=1

    while (log_terms[k] + log(abs(L + (1/(1 - L*(L*L)) + 1/(expm1(log_terms[k] - log_terms[k-1]))))) >= log(mpf(2)) + leps and k < M+initial_k):
        k+=1
        log_terms[k] = f(theta, k)
    
    log_sum = logsumexp(log_terms[initial_k:(k+1)])

    t0 = expm1(log_terms[k] - log_terms[k-1])
    return (k-initial_k, logsumexp([log_sum, log(mpf("0.5") + mpf("-0.5")/t0) + (log_terms[k] + L/(1 + (-1)/t0))]))
