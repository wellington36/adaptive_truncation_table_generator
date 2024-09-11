from mpmath import mp, mpf, exp, log, expm1, fabs, loggamma
from utils import *

### utils
def logsumexp(a):
    with mp.extradps(10):
        return mp.log(mp.fsum([mp.exp(ai) for ai in a]))

### functions
def log_Z_brute(log_lambda: float, nu: float, M: int):
    k = 2
    log_Z_terms = [mpf(0)] * M

    if (nu == 1):
        print("nu = 1")
    if (nu <= 0):
        print("nu must be positive")
    if (nu == mpf("inf")):
        print("nu must be finite")
    
    log_Z_terms[1] = log_lambda
    while (k < M-1):
        log_Z_terms[k] = (k-1) * log_lambda - nu * loggamma(k)
        k+=1
        if k % 10000 == 0:
            print(f"{(k)/10000}%")
    
    log_Z = logsumexp(log_Z_terms)

    return log_Z

mp.dps = 400

mu = mpf(1000)
nu = mpf(0.001)
lamb = mu**nu

print(log_Z_brute(log(lamb), nu, 10**6))


