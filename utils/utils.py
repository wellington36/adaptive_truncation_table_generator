from mpmath import mp, exp, log, fsum, fmul


def logsumexp(a: list):
    with mp.extradps(100):  # Increase precision
        a_max = max(a)

        return a_max + log(fsum([exp(ai - a_max) for ai in a]))

def fma(a, b, c):
    with mp.extradps(100):  # Extra precision for more accurate result
        return fmul(a, b) + c

def logdiffexp(a, b):
    if b > a:
        a, b = b, a

    result = a + log(1 - exp(b - a))

    return result