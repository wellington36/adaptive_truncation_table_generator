from src.bounding_pairs_mp import bounding_pairs_mp
from src.sequential_mp import sequential_mp
from src.fixed_mp import fixed_mp

# pip install /path/to/local/clone
import pybind_stan_fns as psf

from rpy2.robjects.packages import importr
from mpmath import mp, mpf, log, exp, loggamma
from utils.utils import logdiffexp
from tabulate import tabulate
import math


# function of each term in log scale
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
        return (mpf(k)-1) * theta[0] - theta[1] * loggamma(mpf(k))


if __name__ == "__main__":
    mp.dps = 200

    mu = [mpf(10), mpf(100), mpf(1000), mpf(10000)]
    nu = [mpf("0.1"), mpf("0.01"), mpf("0.001"), mpf("0.0001")]
    lamb = [mu[i]**nu[i] for i in range(0,4)]
    loglamb = [log(x) for x in lamb]
    initial_k = 1
    M = [10**4, 10**5, 10**5, 3*10**5]

    # error = 2.2x10^-10
    error = mpf(2)**mpf(-52) * 10**6
    error_minus_10 = []
    for i in range(len(mu)):
        fixed_value = fixed_mp(f, (loglamb[i], nu[i]), M[i], initial_k=initial_k)[1]
        bp_iter = bounding_pairs_mp(f, (loglamb[i], nu[i]), M[i], mpf(0), error, initial_k=initial_k)[1]
        sequential_iter = sequential_mp(f, (loglamb[i],nu[i]), M[i], error, initial_k=initial_k)[1]

        error_minus_10.append([exp(logdiffexp(sequential_iter, fixed_value)), exp(logdiffexp(bp_iter, fixed_value))])

    # error = 2.2x10^-16
    error = mpf(2)**mpf(-52)
    error_minus_16 = []
    for i in range(len(mu)):
        fixed_value = fixed_mp(f, (loglamb[i], nu[i]), M[i], initial_k=initial_k)[1]
        bp_iter = bounding_pairs_mp(f, (loglamb[i], nu[i]), M[i], mpf(0), error, initial_k=initial_k)[1]
        sequential_iter = sequential_mp(f, (loglamb[i],nu[i]), M[i], error, initial_k=initial_k)[1]

        error_minus_16.append([exp(logdiffexp(sequential_iter, fixed_value)), exp(logdiffexp(bp_iter, fixed_value))])
    
    # Libraries
    comp_reg = importr('COMPoissonReg')
    brms_fixed_comp = psf.expose('stan/comp_Z_brms_fixed.stan')

    def dcmp_in_log_scale(x, lambda_, nu):
        # Call the dcmp function with log=TRUE in R
        result = comp_reg.dcmp(x, lambda_, nu)
        return float(result[0])  # Return the first result as a Python float

    libraries = []
    for i in range(len(mu)):
        fixed_value = fixed_mp(f, (loglamb[i], nu[i]), M[i], initial_k=initial_k)[1]

        brms = brms_fixed_comp.log_Z_com_poisson(float(math.log(mu[i])), float(nu[i]))
        dcmp = log(list(comp_reg.dcmp(0, float(lamb[i]), float(nu[i])))[0])

        libraries.append([exp(logdiffexp(fixed_value, brms)), exp(logdiffexp(fixed_value, -1*mpf(dcmp)))])

    
    # Organize in a table
    data = []
    for idx, (a, b, c) in enumerate(zip(error_minus_10, error_minus_16, libraries), start=1):
        data.append([f"mu={10**idx} | nu={10**-idx}", a[0], a[1], b[0], b[1], c[0], c[1]])

    headers = ["", "2.2x10^-10|Sequential", "2.2x10^-10|BP", "2.2x10^-16|Sequential", "2.2x10^-16|BP", "brms", "COMPoissonReg"]

    print(tabulate(data, headers, tablefmt="fancy_grid"))