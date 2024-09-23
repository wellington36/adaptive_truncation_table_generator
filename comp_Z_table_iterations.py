from src.bounding_pairs_mp import bounding_pairs_mp
from src.sequential_mp import sequential_mp
from tabulate import tabulate
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

    mu = [mpf(10), mpf(100), mpf(1000), mpf(10000)]
    nu = [mpf("0.1"), mpf("0.01"), mpf("0.001"), mpf("0.0001")]
    lamb = [mu[i]**nu[i] for i in range(0,4)]
    loglamb = [log(x) for x in lamb]
    M = [10**4, 10**5, 10**5, 3*10**5]

    # 2.2x10^-10
    error = mpf(2)**mpf(-52) * 10**6
    error_minus_10 = []
    for i in range(len(mu)):
        bp_iter = bounding_pairs_mp(f, (loglamb[i], nu[i]), M[i], mpf(0), error, initial_k=1)[0]
        sequential_iter = sequential_mp(f, (loglamb[i],nu[i]), M[i], error, initial_k=1)[0]

        error_minus_10.append([sequential_iter, bp_iter])

    # 2.2x10^-16
    error = mpf(2)**mpf(-52)
    error_minus_16 = []
    for i in range(len(mu)):
        bp_iter = bounding_pairs_mp(f, (loglamb[i], nu[i]), M[i], mpf(0), error, initial_k=1)[0]
        sequential_iter = sequential_mp(f, (loglamb[i],nu[i]), M[i], error, initial_k=1)[0]

        error_minus_16.append([sequential_iter, bp_iter])
    
    # Organize in a table
    data = []
    for idx, (a, b) in enumerate(zip(error_minus_10, error_minus_16), start=1):
        data.append([f"mu={10**idx} | nu={10**-idx}", a[0], a[1], b[0], b[1]])

    headers = ["", "2.2x10^-10|Sequential", "2.2x10^-10|BP", "2.2x10^-16|Sequential", "2.2x10^-16|BP"]

    print(tabulate(data, headers, tablefmt="fancy_grid"))