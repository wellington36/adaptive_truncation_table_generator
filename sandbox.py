from mpmath import mp, mpf, exp, log, expm1, fabs
from utils.utils import *
from math import lgamma
import mpmath


### functions
def log_Z_brute(log_lambda: float, nu: float, M: int):
    k = 2
    log_Z_terms = [log(mpf(0))] * (M+1)

    if (nu == 1):
        print("nu = 1")
    if (nu <= 0):
        print("nu must be positive")
    if (nu == mpf("inf")):
        print("nu must be finite")
    
    log_Z_terms[0] = log(mpf(0))
    log_Z_terms[1] = mpf(0)
    log_Z_terms[2] = log_lambda
    while (k < M):
        k+=1
        log_Z_terms[k] = (k-1) * log_lambda - nu * lgamma(k)

    
    log_Z = logsumexp(log_Z_terms)

    return log_Z

def log_Z_value_bp(log_lambda: float, nu: float, M: int, eps: float):
    k = 2
    leps = log(eps)
    log_Z_terms = [log(mpf(0))] * M

    if (nu == 1):
        print("nu = 1")
    if (nu <= 0):
        print("nu must be positive")
    if (nu == mpf("inf")):
        print("nu must be finite")
    
    log_Z_terms[0] = log(mpf(0))
    log_Z_terms[1] = 0
    log_Z_terms[2] = log_lambda
    while (log_Z_terms[k] + log(abs(1 + 1/expm1(log_Z_terms[k] - log_Z_terms[k-1]))) >= log(mpf(2)) + leps and k < M):
        k+=1
        log_Z_terms[k] = (k-1) * log_lambda - nu * lgamma(k)
    
    log_Z = logsumexp(log_Z_terms)

    return (k-1, logsumexp([log_Z, log_Z_terms[k] + log(0.5 - 0.5/expm1(log_Z_terms[k] - log_Z_terms[k-1]))]))

mp.dps = 200

mu = [mpf(10), mpf(100), mpf(1000), mpf(10000)]
nu = [mpf("0.1"), mpf("0.01"), mpf("0.001"), mpf("0.0001")]
lamb = [mu[i]**nu[i] for i in range(0,4)]
M = 10**5

for j in [8, 7, 6, 5]:
    eps = mpf(2**(-52)) * 10**j

    for i in range(0, 4):
        brute_Z = log_Z_brute(log(lamb[i]), nu[i], M)
        (bounding_k, bounding_Z) = log_Z_value_bp(log(lamb[i]), nu[i], M, eps)
    
        error = mpmath.nstr(abs(exp(logdiffexp(brute_Z, bounding_Z))), n=5, strip_zeros=False)
        str_eps = mpmath.nstr(eps, n=2, strip_zeros=False)

        print(f"eps = {str_eps} | mu = {mu[i]}      | nu = {nu[i]}      | error = {error}       | k = {bounding_k}")



