from src.bounding_pairs_mp import bounding_pairs_mp
from src.sum_to_threshold_mp import sum_to_threshold_mp
from src.fixed_mp import fixed_mp

from mpmath import mpf, log, exp, loggamma
from utils.utils import logdiffexp
from tabulate import tabulate


# function of each term in log scale
def f(theta: tuple, k: int):
    lamb = mpf(theta[0])
    r = mpf(theta[1])
    k = mpf(k)

    return k * log(lamb) - lamb - loggamma(k-r+1)

def f_true(theta):
    return theta[1] * log(theta[0])

def ratio_bounding(lamb, r, eps, M, brute):
    check1 = 0
    check2 = 0
    check3 = 0

    k, aprox = bounding_pairs_mp(f, (lamb, mpf(r)), L=0, eps=eps, M=M, initial_k=r)

    print(k)
    print(exp(logdiffexp(aprox, f_true((lamb, mpf(r))))))

    if exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= eps:
        check1 = 1
    
    if exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= exp(logdiffexp(f_true((lamb, mpf(r))), brute)):
        check2 = 1
        
    if (exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= eps) or (exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= exp(logdiffexp(f_true((lamb, mpf(r))), brute))):
        check3 = 1

    return (check1, check2, check3)

def ratio_threshold(lamb, r, eps, M, brute):
    check1 = 0
    check2 = 0
    check3 = 0

    _, aprox = sum_to_threshold_mp(f, (lamb, mpf(r)), L=0, eps=eps, M=M, initial_k=r)

    if exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= eps:
        check1 = 1
    
    if exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= exp(logdiffexp(f_true((lamb, mpf(r))), brute)):
        check2 = 1
        
    if (exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= eps) or (exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= exp(logdiffexp(f_true((lamb, mpf(r))), brute))):
        check3 = 1

    return (check1, check2, check3)

def ratio_fixed(lamb, r, eps, m, M, brute):
    check1 = 0
    check2 = 0
    check3 = 0

    if M == m:
        return (exp(logdiffexp(brute, f_true((lamb, mpf(r)))) <= eps), 1, 1)

    _, aprox = fixed_mp(f, (lamb, mpf(r)), M=m, initial_k=r)

    if exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= eps:
        check1 = 1
    
    if exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= exp(logdiffexp(f_true((lamb, mpf(r))), brute)):
        check2 = 1
        
    if (exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= eps) or (exp(logdiffexp(aprox, f_true((lamb, mpf(r))))) <= exp(logdiffexp(f_true((lamb, mpf(r))), brute))):
        check3 = 1

    return (check1, check2, check3)


if __name__ == "__main__":
    #mp.dps = 100
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
            _, brute = fixed_mp(f, (lamb, mpf(r)), M=M, initial_k=r)

            for eps in epslist:
                print(f"Evaluating: lambda={lamb}, r={r}, eps={float(eps)}...")
                
                ratios_bounding = [sum(i) for i in zip(ratios_bounding, ratio_bounding(lamb, r, eps, M, brute))]
                print(ratios_bounding)
                ratios_threshold = [sum(i) for i in zip(ratios_threshold, ratio_threshold(lamb, r, eps, M, brute))]
                print(ratios_threshold)
                ratios_cap_m = [sum(i) for i in zip(ratios_cap_m, ratio_fixed(lamb, r, eps, m, M, brute))]
                print(ratios_cap_m)
                ratios_cap_M = [sum(i) for i in zip(ratios_cap_M, ratio_fixed(lamb, r, eps, M, M, brute))]
                print(ratios_cap_M)

    Steps = len(lamblist) * len(rlist) * len(epslist)

    ratios_bounding = [x/Steps for x in ratios_bounding]
    ratios_threshold = [x/Steps for x in ratios_threshold]
    ratios_cap_m = [x/Steps for x in ratios_cap_m]
    ratios_cap_M = [x/Steps for x in ratios_cap_M]
    methods = ["Error-bounding pairs", "Threshold", "Cap = 1 x 10^3", "Cap = 5 x 10^5"]

    ratios = [ratios_bounding, ratios_threshold, ratios_cap_m, ratios_cap_M]

    data = []
    for idx, (a) in enumerate(ratios):
        data.append([f"{methods[idx]}", a[0], a[1], a[2]])

    headers = ["", "Error", "Error with M", "Either"]

    print(tabulate(data, headers, tablefmt="fancy_grid"))
