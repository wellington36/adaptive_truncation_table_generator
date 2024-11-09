from src.bounding_pairs_mp import bounding_pairs_mp
from src.sum_to_threshold_mp import sum_to_threshold_mp
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

    theta   : (lambda, r)
    k       : k-th term
    """
    lamb = mpf(theta[0])
    r = mpf(theta[1])
    k = mpf(k)

    return k * log(lamb) - lamb - loggamma(k-r+1)

def f_true(theta):
    return theta[1] * log(theta[0])

def ratio_bounding(lamb, r, eps, M):
    check1 = 0
    check2 = 0
    check3 = 0

    _, aprox = bounding_pairs_mp(f, (lamb, mpf(r)), L=0, eps=eps, M=M, initial_k=r)
    _, brute = fixed_mp(f, (lamb, mpf(r)), M=M, initial_k=r)

    if exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= eps:
        check1 = 1
    
    if exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= exp(logdiffexp(f_true((lamb, mpf(r))), brute)):
        check2 = 1
        
    if (exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= eps) or (exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= exp(logdiffexp(f_true((lamb, mpf(r))), brute))):
        check3 = 1

    return (check1, check2, check3)

def ratio_threshold(lamb, r, eps, M):
    check1 = 0
    check2 = 0
    check3 = 0

    _, aprox = sum_to_threshold_mp(f, (lamb, mpf(r)), L=0, eps=eps, M=M, initial_k=r)
    _, brute = fixed_mp(f, (lamb, mpf(r)), M=M, initial_k=r)

    if exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= eps:
        check1 = 1
    
    if exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= exp(logdiffexp(f_true((lamb, mpf(r))), brute)):
        check2 = 1
        
    if (exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= eps) or (exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= exp(logdiffexp(f_true((lamb, mpf(r))), brute))):
        check3 = 1

    return (check1, check2, check3)

def ratio_fixed(lamb, r, eps, m, M):
    check1 = 0
    check2 = 0
    check3 = 0

    if M == m:
        return (logdiffexp(fixed_mp(f, (lamb, mpf(r)), M=M, initial_k=r)[1], f_true((lamb, mpf(r)))) <= eps, 1, 1)

    _, aprox = fixed_mp(f, (lamb, mpf(r)), M=m, initial_k=r)
    _, brute = fixed_mp(f, (lamb, mpf(r)), M=M, initial_k=r)

    if exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= eps:
        check1 = 1
    
    if exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= exp(logdiffexp(f_true((lamb, mpf(r))), brute)):
        check2 = 1
        
    if (exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= eps) or (exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= exp(logdiffexp(f_true((lamb, mpf(r))), brute))):
        check3 = 1

    return (check1, check2, check3)

if __name__ == "__main__":
    mpf.dps = 200
    machine_eps = mpf(2)**mpf(-52)

    lamblist = [mpf("0.5"), mpf(1), mpf(10), mpf(1000)]
    rlist = [2, 5, 10]
    epslist = [machine_eps, machine_eps*10, machine_eps*10**4]
    m = 1*10**3
    M = 5*10**5

    ratios_bounding = [0, 0, 0]
    ratios_threshold = [0, 0, 0]
    ratios_cap_m = [0, 0, 0]
    ratios_cap_M = [0, 0, 0]

    for lamb in lamblist:
        for r in rlist:
            for eps in epslist:
                print(f"Evaluating: lambda={lamb}, r={r}, eps={float(eps)}...")
                
                ratios_bounding = [sum(i) for i in zip(ratios_bounding, ratio_bounding(lamb, r, eps, M))]
                print(ratios_bounding)
                ratios_threshold = [sum(i) for i in zip(ratios_threshold, ratio_threshold(lamb, r, eps, M))]
                print(ratios_threshold)
                ratios_cap_m = [sum(i) for i in zip(ratios_cap_m, ratio_fixed(lamb, r, eps, m, M))]
                print(ratios_cap_m)
                ratios_cap_M = [sum(i) for i in zip(ratios_cap_M, ratio_fixed(lamb, r, eps, M, M))]
                print(ratios_cap_M)

    Steps = len(lamblist) * len(rlist) * len(epslist)

    ratios_bounding = [x/Steps for x in ratios_threshold]
    ratios_threshold = [x/Steps for x in ratios_threshold]
    ratios_cap_m = [x/Steps for x in ratios_cap_m]
    ratios_cap_M = [x/Steps for x in ratios_cap_M]
    methods = ["Error-bounding pairs", "Threshold", "Cap = 1 x 10^3", "Cap = 5 x 10^5"]

    ratios = [ratios_bounding, ratios_threshold, ratios_cap_m, ratios_cap_M]

    data = []
    for idx, (a) in enumerate(zip(ratios)):
        data.append([f"{methods[idx]}", a[0], a[1], a[2]])

    headers = ["", "Error", "Error with M", "Either"]

    print(tabulate(data, headers, tablefmt="fancy_grid"))


    #print(ratio_bounding(lamb[3], r[2], eps[0], M, N))
