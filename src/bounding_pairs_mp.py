from mpmath import mp, mpf, exp, log, expm1
from utils.utils import logsumexp
from math import lgamma
import mpmath

def bounding_pairs_mp(f, theta, M, L, eps, initial_k):
    k = initial_k
    leps = log(eps)
    log_terms = [log(mpf(1))] * M

    log_terms[k] = f(theta, k)
    log_terms[k+1] = f(theta, k+1)
    k+=1

    while (log_terms[k] + log(abs(L + (1/(1 - L*(L*L)) + 1/(expm1(log_terms[k] - log_terms[k-1]))))) >= log(mpf(2)) + leps and k < M-1):
        k+=1
        print(log_terms[k] + log(abs(L + (1/(1 - L*(L*L)) + 1/(expm1(log_terms[k] - log_terms[k-1]))))))
        print(log(mpf(2)) + leps)
        log_terms[k] = f(theta, k)
    
    log_Z = logsumexp(log_terms)

    t0 = expm1(log_terms[k] - log_terms[k-1])
    return (k-initial_k, logsumexp([log_Z, log(mpf("0.5") + mpf("-0.5")/t0) + (log_terms[k] + L/(1 + (-1)/t0))]))
