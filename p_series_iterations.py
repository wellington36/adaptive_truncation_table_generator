from src.integration_bound_mp import integration_bound_mp
from src.sequential_mp import sequential_mp
from mpmath import mpf, log, exp, mp

def f (p, n):
    return log(mpf(1)) - (p * log(mpf(n)))

def g (p, n):
    return (1 - p) * log(mpf(n)) - log(mpf(p) - 1)

mp.dps = 200
pl = [mpf(3), mpf(2), mpf('1.5'), mpf('1.4'), mpf('1.3')]
epsl = [mpf(2)**(mpf(-52))*mpf(10)**mpf(14), mpf(2)**(mpf(-52))*mpf(10)**mpf(12)]
M = [10**5, 10**5, 2*10**6, 5*10**6, 10**7]


for eps in epsl:
    for p in pl:
        integration = integration_bound_mp(f, g, p, M[pl.index(p)], eps, 1)
        sequential = sequential_mp(f, p, M[pl.index(p)], eps, 1)

        print(f"eps = {eps} | p = {p} | integration = {integration[0]} | sequential = {sequential[0]}")
