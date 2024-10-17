from src.bounding_pairs_mp import bounding_pairs_mp
from src.fixed_mp import fixed_mp
from utils.utils import logdiffexp

from mpmath import mp, mpf, log, pi, exp, loggamma, fsum
from tabulate import tabulate
import time


def f(theta: tuple, k: int):
    """
    terms of the normalization contant series

    theta   : (log_lambda, nu)
    k       : k-th term
    """
    return (mpf(k)) * theta[0] - theta[1] * loggamma(mpf(k)+1)

def log_Z_comp_approx_max(log_lambda, nu):
    log_lambda = mpf(log_lambda)
    nu = mpf(nu)

    nu2 = nu**2
    log_common = log(nu) + log_lambda/nu
    resids = [log(mpf(0))] * 4
    lcte = (nu * exp(log_lambda/nu)) - ((nu - 1)/(2*nu)*log_lambda + (nu - 1)/2*log(2*pi) + mpf('0.5')*log(nu))

    c1 = (nu2 - 1) / 24
    c2 = (nu2 - 1) / 1152 * (nu2 + 23)
    c3 = (nu2 - 1) / 414720 * (5*nu2**2 - 298*nu2 + 11237)

    resids[0] = 1
    resids[1] = c1 * exp(-1 * log_common)
    resids[2] = c2 * exp(-2 * log_common)
    resids[3] = c3 * exp(-3 * log_common)

    return lcte + log(fsum(resids))



if __name__ == '__main__':
    mp.dps = 300
    M = 10**5
    L = 0
    eps = mpf(2)**mpf(-52)
    initial_k = 0

    lamb = [mpf('1.5'), mpf('1.5')+eps*10**10, mpf('1.5'), mpf(5), mpf('1.5'), mpf(5)]
    nu = [mpf(1), mpf(1), mpf(1)+eps*10**10, mpf('0.8'), mpf(2), mpf(2)]

    bounding_pairs_error = []
    bounding_pairs_time = []
    approx_error = []
    approx_time = []

    for i in range(len(lamb)):
        print(f"Init: lambda = {float(lamb[i])} | nu = {float(nu[i])}")
        # test bounding pairs
        t0_bound = time.time()
        k, log_Z_bound = bounding_pairs_mp(f, (log(lamb[i]), nu[i]), M, L, eps, initial_k=initial_k)
        t1_bound = time.time() - t0_bound

        t0_bound = time.time()
        _, _ = bounding_pairs_mp(f, (log(lamb[i]), nu[i]), M, L, eps, initial_k=initial_k)
        t2_bound = time.time() - t0_bound

        t0_bound = time.time()
        _, _ = bounding_pairs_mp(f, (log(lamb[i]), nu[i]), M, L, eps, initial_k=initial_k)
        t3_bound = time.time() - t0_bound

        # test approximations
        t0_approx = time.time()
        log_Z_approx = log_Z_comp_approx_max(log(lamb[i]), nu[i])
        t1_approx = time.time() - t0_approx

        t0_approx = time.time()
        _ = log_Z_comp_approx_max(log(lamb[i]), nu[i])
        t2_approx = time.time() - t0_approx

        t0_approx = time.time()
        _ = log_Z_comp_approx_max(log(lamb[i]), nu[i])
        t3_approx = time.time() - t0_approx

        # true
        print("Fixed evaluation...")
        _, fixed = fixed_mp(f, (log(lamb[i]), nu[i]), M, initial_k=initial_k)

        # storage variables
        bounding_pairs_error.append((f"{float(exp(logdiffexp(log_Z_bound, fixed)))} / {k}"))
        approx_error.append(float(exp(logdiffexp(log_Z_approx, fixed))))

        bounding_pairs_time.append((t1_bound + t2_bound + t3_bound)/3)
        approx_time.append((t1_approx + t2_approx + t3_approx)/3)


    # Organize in a table
    data = []
    for idx, (a, b, c, d) in enumerate(zip(bounding_pairs_error, approx_error, bounding_pairs_time, approx_time)):
        data.append([f"lambda={float(lamb[idx])} | nu={float(nu[idx])}", a, b, c*1000, d*1000])

    headers = ["", "Error|Bounding / k", "Error|Approx", "Time(ms)|Bounding", "Time(ms)|Approx"]

    print(tabulate(data, headers, tablefmt="fancy_grid"))
