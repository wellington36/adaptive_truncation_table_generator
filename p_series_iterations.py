from src.integration_bound_mp import integration_bound_mp
from src.sequential_mp import sequential_mp
from mpmath import mpf, log, exp, mp
from tabulate import tabulate

def f (p, n):
    return log(mpf(1)) - (p * log(mpf(n)))

def g (p, n):
    return (1 - p) * log(mpf(n)) - log(mpf(p) - 1)

mp.dps = 200
pl = [mpf(3), mpf(2), mpf('1.6'), mpf('1.5'), mpf('1.4')]
epsl = [mpf(2)**(mpf(-52))*mpf(10)**mpf(14), mpf(2)**(mpf(-52))*mpf(10)**mpf(12)]
M = [10**5, 10**5, 2*10**6, 5*10**6, 5*10**6]
initial_k = 1


terms_error_minus_2 = []
for p in pl:
    print(f"Evaluating: eps = {float(mpf(2)**(mpf(-52))*mpf(10)**mpf(14))} and p = {float(p)}")

    eps = mpf(2)**(mpf(-52))*mpf(10)**mpf(14)
    _, integration = integration_bound_mp(f, g, p, M[pl.index(p)], eps, initial_k)
    _, sequential = sequential_mp(f, p, M[pl.index(p)], eps, initial_k)

    terms_error_minus_2.append((sequential, integration))

terms_error_minus_4 = []
for p in pl:
    print(f"Evaluating: eps = {float(mpf(2)**(mpf(-52))*mpf(10)**mpf(12))} and p = {float(p)}")

    eps = mpf(2)**(mpf(-52))*mpf(10)**mpf(12)
    _, integration = integration_bound_mp(f, g, p, M[pl.index(p)], eps, initial_k)
    _, sequential = sequential_mp(f, p, M[pl.index(p)], eps, initial_k)

    terms_error_minus_4.append((sequential, integration))

    data = []
    for idx, (a, b) in enumerate(zip(terms_error_minus_2, terms_error_minus_4)):
        data.append([f"p={pl[idx]}", a[0], a[1], b[0], b[1]])

    headers = ["", "2.2x10^-2|Sequential", "2.2x10^-2|Integration", "2.2x10^-4|Sequential", "2.2x10^-4|Integration"]

    print(tabulate(data, headers, tablefmt="fancy_grid"))