from src.bounding_pairs_mp import bounding_pairs_mp
from src.sequential_mp import sequential_mp
from src.brute_mp import brute_mp
from tabulate import tabulate
from mpmath import mp, mpf, log, exp
from math import lgamma
from rpy2.robjects.packages import importr
from utils.utils import logdiffexp


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
    #M = [10**4, 10**4, 10**5, 10**5]

    # 2.2x10^-10
    error = mpf(2)**mpf(-52) * 10**6
    error_minus_10 = []
    for i in range(len(mu)):
        brute_value = brute_mp(f, (loglamb[i], nu[i]), M[i], initial_k=1)[1]
        bp_iter = bounding_pairs_mp(f, (loglamb[i], nu[i]), M[i], mpf(0), error, initial_k=1)[1]
        sequential_iter = sequential_mp(f, (loglamb[i],nu[i]), M[i], error, initial_k=1)[1]

        error_minus_10.append([exp(logdiffexp(sequential_iter, brute_value)), exp(logdiffexp(bp_iter, brute_value))])

    # 2.2x10^-16
    error = mpf(2)**mpf(-52)
    error_minus_16 = []
    for i in range(len(mu)):
        brute_value = brute_mp(f, (loglamb[i], nu[i]), M[i], initial_k=1)[1]
        bp_iter = bounding_pairs_mp(f, (loglamb[i], nu[i]), M[i], mpf(0), error, initial_k=1)[1]
        sequential_iter = sequential_mp(f, (loglamb[i],nu[i]), M[i], error, initial_k=1)[1]

        error_minus_16.append([exp(logdiffexp(sequential_iter, brute_value)), exp(logdiffexp(bp_iter, brute_value))])
    
    # Libraries
    comp_reg = importr('COMPoissonReg')

    def dcmp_in_log_scale(x, lambda_, nu):
        # Call the dcmp function with log=TRUE in R
        result = comp_reg.dcmp(x, lambda_, nu)
        return float(result[0])  # Return the first result as a Python float

    libraries = []
    for i in range(len(mu)):
        brute_value = brute_mp(f, (loglamb[i], nu[i]), M[i], initial_k=1)[1]
        dcmp = log(list(comp_reg.dcmp(0, float(lamb[i]), float(nu[i])))[0])

        libraries.append([abs(exp(logdiffexp(brute_value, -1*mpf(dcmp))))])

    
    # Organize in a table
    data = []
    for idx, (a, b, c) in enumerate(zip(error_minus_10, error_minus_16, libraries), start=1):
        data.append([f"mu={10**idx} | nu={10**-idx}", a[0], a[1], b[0], b[1], c[0]])

    headers = ["", "2.2x10^-10|Sequential", "2.2x10^-10|BP", "2.2x10^-16|Sequential", "2.2x10^-16|BP", "COMPoissonReg"]

    print(tabulate(data, headers, tablefmt="fancy_grid"))