from src.bounding_pairs_mp import bounding_pairs_mp
from mpmath import mp, mpf, log
from math import lgamma

def f(theta: tuple, k: int):
    """
    terms of the normalization contant series

    theta   : (log_lambda, nu)
    k       : k-th term
    """

    if k == 1:
        return mpf(0)
    elif (k == 2):
        return theta[0]
    else:
        return (k-1) * theta[0] - theta[1] * lgamma(k)


if __name__ == "__main__":
    mp.dps = 200
    mu = mpf(10)
    nu = mpf('0.1')
    lamb = mu**nu

    theta = (log(lamb), nu)
    
    print("Hello, world!")
    print(bounding_pairs_mp(f, theta, 10**4, mpf(0), eps=mpf(2)**mpf(-52), initial_k=1))