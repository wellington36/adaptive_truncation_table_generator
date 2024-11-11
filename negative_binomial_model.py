from src.bounding_pairs_mp import bounding_pairs_mp
from src.sum_to_threshold_mp import sum_to_threshold_mp
from src.fixed_mp import fixed_mp

from mpmath import mp, mpf, log, exp, binomial
from utils.utils import logdiffexp
from tabulate import tabulate


def f(theta: tuple, y: int):
    phi, mu, x, eta = theta

    term1 = log(binomial(y + phi - 1, y))     # Log of binomial coefficient
    term2 = y * log(mu / (mu + phi))          # Log of (mu / (mu + phi))^y
    term3 = phi * log(phi / (mu + phi))       # Log of (phi / (mu + phi))^phi
    term4 = log(binomial(y, x))               # Log of binomial coefficient (y choose x)
    term5 = x * log(eta)                      # Log of eta^x
    term6 = (y - x) * log(1 - eta)            # Log of (1 - eta)^(y - x)

    return term1 + term2 + term3 + term4 + term5 + term6

def f_true(theta):
    phi, mu, x, eta = theta

    term1 = log(binomial(x + phi - 1, x))
    term2 = x * log((eta * mu)/(eta * mu + phi))
    term3 = phi * log(phi/(eta * mu + phi))

    return term1 + term2 + term3


def ratio_bounding(phi, mu, x, eta, eps, M, L, brute):
    check1 = 0
    check2 = 0
    check3 = 0

    k, aprox = bounding_pairs_mp(f, (phi, mu, mpf(x), eta), L=L, eps=eps, M=M, initial_k=x)

    print(k)
    print(exp(logdiffexp(aprox, f_true((phi, mu, mpf(x), eta)))))

    if exp(logdiffexp(aprox, f_true((phi, mu, mpf(x), eta)))) <= eps:
        check1 = 1
    
    if exp(logdiffexp(aprox, f_true((phi, mu, mpf(x), eta)))) <= exp(logdiffexp(f_true((phi, mu, mpf(x), eta)), brute)):
        check2 = 1
        
    if check1 == 1 or check2 == 1:
        check3 = 1

    return (check1, check2, check3)

def ratio_threshold(phi, mu, x, eta, eps, M, L, brute):
    check1 = 0
    check2 = 0
    check3 = 0

    k, aprox = sum_to_threshold_mp(f, (phi, mu, mpf(x), eta), L=L, eps=eps, M=M, initial_k=x)

    if exp(logdiffexp(aprox, f_true((phi, mu, mpf(x), eta)))) <= eps:
        check1 = 1
    
    if exp(logdiffexp(aprox, f_true((phi, mu, mpf(x), eta)))) <= exp(logdiffexp(f_true((phi, mu, mpf(x), eta)), brute)):
        check2 = 1
        
    if check1 == 1 or check2 == 1:
        check3 = 1

    return (check1, check2, check3)

def ratio_fixed(phi, mu, x, eta, eps, m, M, brute):
    check1 = 0
    check2 = 0
    check3 = 0

    if M == m:
        return (exp(logdiffexp(brute, f_true((phi, mu, mpf(x), eta)))) <= eps, 1, 1)

    _, aprox = fixed_mp(f, (phi, mu, mpf(x), eta), M=m, initial_k=x)

    if exp(logdiffexp(aprox, f_true((phi, mu, mpf(x), eta)))) <= eps:
        check1 = 1
    
    if exp(logdiffexp(aprox, f_true((phi, mu, mpf(x), eta)))) <= exp(logdiffexp(f_true((phi, mu, mpf(x), eta)), brute)):
        check2 = 1
        
    if check1 or check2:
        check3 = 1

    return (check1, check2, check3)


if __name__ == "__main__":
    #mp.dps = 100
    machine_eps = mpf(2)**mpf(-52)

    phil = [mpf('0.1'), mpf('0.5'), mpf(1), mpf(10)]
    mul = [mpf(1), mpf(10), mpf(100)]
    xl = [0, 5, 10]
    etal = [mpf('0.01'), mpf('0.1'), mpf('0.5'), mpf('0.75')]
    epsl = [machine_eps, machine_eps*10, machine_eps*10**4]
    m = 1*10**3
    M = 5*10**5

    ratios_bounding_less = [0, 0, 0]
    ratios_threshold_less = [0, 0, 0]
    ratios_cap_m_less = [0, 0, 0]
    ratios_cap_M_less = [0, 0, 0]
    ratios_bounding_great = [0, 0, 0]
    ratios_threshold_great = [0, 0, 0]
    ratios_cap_m_great = [0, 0, 0]
    ratios_cap_M_great = [0, 0, 0]

    count_of_less = 0
    count_of_great = 0

    for phi in phil:
        for mu in mul:
            for x in xl:
                for eta in etal:
                    print(f"Evaluating: theta=({phi, mu, x, eta})...")
                    
                    _, brute = fixed_mp(f, (phi, mu, x, eta), M=M, initial_k=x)
                    L = (mu/(mu + phi)) * (1 - eta)
                    print(f"True: {exp(f_true((phi, mu, mpf(x), eta)))}")

                    for eps in epsl:
                        if L < 1/2:
                            count_of_less+=1

                            ratios_bounding_less = [sum(i) for i in zip(ratios_bounding_less, ratio_bounding(phi, mu, x, eta, eps, M, L, brute))]
                            print(ratios_bounding_less)
                            ratios_threshold_less = [sum(i) for i in zip(ratios_threshold_less, ratio_threshold(phi, mu, x, eta, eps, M, L, brute))]
                            print(ratios_threshold_less)
                            ratios_cap_m_less = [sum(i) for i in zip(ratios_cap_m_less, ratio_fixed(phi, mu, x, eta, eps, m, M, brute))]
                            print(ratios_cap_m_less)
                            ratios_cap_M_less = [sum(i) for i in zip(ratios_cap_M_less, ratio_fixed(phi, mu, x, eta, eps, M, M, brute))]
                            print(ratios_cap_M_less)
                        else:
                            count_of_great+=1

                            ratios_bounding_great = [sum(i) for i in zip(ratios_bounding_great, ratio_bounding(phi, mu, x, eta, eps, M, L, brute))]
                            print(ratios_bounding_great)
                            ratios_threshold_great = [sum(i) for i in zip(ratios_threshold_great, ratio_threshold(phi, mu, x, eta, eps, M, L, brute))]
                            print(ratios_threshold_great)
                            ratios_cap_m_great = [sum(i) for i in zip(ratios_cap_m_great, ratio_fixed(phi, mu, x, eta, eps, m, M, brute))]
                            print(ratios_cap_m_great)
                            ratios_cap_M_great = [sum(i) for i in zip(ratios_cap_M_great, ratio_fixed(phi, mu, x, eta, eps, M, M, brute))]
                            print(ratios_cap_M_great)


    ratios_bounding_less = [x/count_of_less for x in ratios_bounding_less]
    ratios_threshold_less = [x/count_of_less for x in ratios_threshold_less]
    ratios_cap_m_less = [x/count_of_less for x in ratios_cap_m_less]
    ratios_cap_M_less = [x/count_of_less for x in ratios_cap_M_less]
    ratios_bounding_great = [x/count_of_great for x in ratios_bounding_great]
    ratios_threshold_great = [x/count_of_great for x in ratios_threshold_great]
    ratios_cap_m_great = [x/count_of_great for x in ratios_cap_m_great]
    ratios_cap_M_great = [x/count_of_great for x in ratios_cap_M_great]
    methods = ["Error-bounding pairs | No", "Threshold | No",
               "Cap = 1 x 10^3 | No", "Cap = 5 x 10^5 | No",
               "Error-bounding pairs | Yes", "Threshold | Yes",
               "Cap = 1 x 10^3 | Yes", "Cap = 5 x 10^5 | Yes"]

    ratios = [ratios_bounding_less, ratios_threshold_less, ratios_cap_m_less, ratios_cap_M_less,
              ratios_bounding_great, ratios_threshold_great, ratios_cap_m_great, ratios_cap_M_great]

    data = []
    for idx, (a) in enumerate(ratios):
        data.append([f"{methods[idx]}", a[0], a[1], a[2]])

    headers = ["Method | L > 1/2", "Error", "Error with M", "Either"]

    print(tabulate(data, headers, tablefmt="fancy_grid"))
