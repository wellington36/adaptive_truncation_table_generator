from src.bounding_pairs_mp import bounding_pairs_mp
from src.fixed_mp import fixed_mp

from mpmath import mp, mpf, log, exp
from utils.utils import logdiffexp
from random import random
from math import lgamma
from time import time


# function of each term in log scale
def f(theta: tuple, k: int):
    """
    terms of the normalization contant series

    theta   : (log_lambda, nu)
    k       : k-th term
    """
    theta = [mpf(x) for x in theta]

    if k == 1:
        return mpf(0)
    elif (k == 2):
        return theta[0]
    else:
        return (mpf(k)-1) * theta[0] - theta[1] * lgamma(mpf(k))

if __name__ == "__main__":
    mp.dps = 400
    M = 10**7
    eps = mpf(2)**mpf(-52)
    initial_k = 1
    L = mpf(0)

    mu_l = 0
    mu_u = 1000
    nu_l = 0
    nu_u = 100

    F = 0
    N = 10000

    for i in range(1, N+1):
        t0 = time()
        mu = mu_u * random() + mu_l
        nu = nu_u * random() + nu_l

        lamb = mu**nu
        loglamb = log(lamb)

        theta = (loglamb, nu)

        k, bp_value = bounding_pairs_mp(f, theta, M, L, eps, initial_k)

        fixed_value = fixed_mp(f, theta, M=k*10, initial_k=initial_k)[1]

        if (i % 50 == 0) and exp(logdiffexp(bp_value, fixed_value)) < eps:
            print(f"Pass test {i}/{N}. (k = {k} | time = {(time()-t0):.2f}s)")
        
        if exp(logdiffexp(bp_value, fixed_value)) > eps:
            print(f"Fail test {i}: theta = {theta} | bp_value = {bp_value} | fixed_value = {fixed_value} | error = {float(exp(logdiffexp(bp_value, fixed_value)))}.")
            F+=1


    print(f"Tests faileds {F/N*100:.2f}% ({F}/{N}).")