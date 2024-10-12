from mpmath import mp, exp, log, fsum


def logsumexp(a: list):
    with mp.extradps(100):  # Increase precision
        a_max = max(a)

        return a_max + log(fsum([mp.exp(ai - a_max) for ai in a]))

def logdiffexp(a, b):
    if b > a:
        a, b = b, a

    result = a + log(1 - exp(b - a))

    return result