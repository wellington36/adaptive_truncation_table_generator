from src.bounding_pairs_mp import bounding_pairs_mp
from src.sequential_mp import sequential_mp
from src.brute_mp import brute_mp

# pip install /path/to/local/clone
import pybind_stan_fns as psf

from rpy2.robjects.packages import importr
from mpmath import mp, mpf, log, exp
from utils.utils import logdiffexp
from tabulate import tabulate
from math import lgamma
import math


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

    mu = [mpf('1.5'), mpf(1), mpf(2), mpf(5), mpf(7), mpf(9), mpf(10)]
    nu = [mpf(1), mpf(1), mpf(2), mpf('0.7'), mpf('0.5'), mpf('0.3'), mpf('0.2')]
    lamb = [mu[i]**nu[i] for i in range(len(mu))]
    loglamb = [log(x) for x in lamb]
    initial_k = 1

    print(lamb)
    print(loglamb)
    M = [10**5]*len(mu)
    #M = [10**4]*len(mu)  # test

    # error = 2.2x10^-10
    error = mpf(2)**mpf(-52) * 10**12
    error_minus_10 = []
    for i in range(len(mu)):
        brute_value = brute_mp(f, (loglamb[i], nu[i]), M[i], initial_k=initial_k)[1]
        k, bp_iter = bounding_pairs_mp(f, (loglamb[i], nu[i]), M[i], mpf(0), error, initial_k=initial_k)
        sequential_iter = sequential_mp(f, (loglamb[i],nu[i]), M[i], error, initial_k=initial_k)[1]

        error_minus_10.append([exp(logdiffexp(sequential_iter, brute_value)), f"{float(exp(logdiffexp(bp_iter, brute_value)))} || {k}"])

    # error = 2.2x10^-16
    error = mpf(2)**mpf(-52) * 10**10
    error_minus_16 = []
    for i in range(len(mu)):
        brute_value = brute_mp(f, (loglamb[i], nu[i]), M[i], initial_k=initial_k)[1]
        k, bp_iter = bounding_pairs_mp(f, (loglamb[i], nu[i]), M[i], mpf(0), error, initial_k=initial_k)
        sequential_iter = sequential_mp(f, (loglamb[i],nu[i]), M[i], error, initial_k=initial_k)[1]

        error_minus_16.append([exp(logdiffexp(sequential_iter, brute_value)), f"{float(exp(logdiffexp(bp_iter, brute_value)))} || {k}"])
    
    # Libraries
    comp_reg = importr('COMPoissonReg')
    brms_fixed_comp = psf.expose('stan/comp_Z_brms_fixed.stan')

    def dcmp_in_log_scale(x, lambda_, nu):
        # Call the dcmp function with log=TRUE in R
        result = comp_reg.dcmp(x, lambda_, nu)
        return float(result[0])  # Return the first result as a Python float

    libraries = []
    for i in range(len(mu)):
        brute_value = brute_mp(f, (loglamb[i], nu[i]), M[i], initial_k=initial_k)[1]

        brms = brms_fixed_comp.log_Z_com_poisson(float(math.log(mu[i])), float(nu[i]))
        dcmp = log(list(comp_reg.dcmp(0, float(lamb[i]), float(nu[i])))[0])

        libraries.append([exp(logdiffexp(brute_value, brms)), exp(logdiffexp(brute_value, -1*mpf(dcmp)))])

    
    # Organize in a table
    data = []
    for idx, (a, b, c) in enumerate(zip(error_minus_10, error_minus_16, libraries)):
        data.append([f"mu={mu[idx]} | nu={nu[idx]}", a[0], a[1], b[0], b[1], c[0], c[1]])

    headers = ["", "2.2x10^-4|Sequential", "2.2x10^-4|BP || k", "2.2x10^-6|Sequential", "2.2x10^-6|BP || k", "brms", "COMPoissonReg"]

    print(tabulate(data, headers, tablefmt="fancy_grid"))