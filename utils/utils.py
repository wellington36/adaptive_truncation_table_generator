from mpmath import mp, mpf, exp, log, expm1, fabs

def logsumexp(a):
    with mp.extradps(100):
        return mp.log(mp.fsum([mp.exp(ai) for ai in a]))

def logdiffexp(a, b):
    if b > a:
        a, b = b, a

    result = a + log(1 - exp(b - a))

    return result