from mpmath import mpf, log, expm1
from utils.utils import logsumexp


def bounding_pairs_mp(f, theta, M, L, eps, initial_k):
    k = initial_k
    leps = log(eps)
    log_terms = [log(mpf(0))] * (M+initial_k)

    if L == 0:
        is_decreasing = True
    else:
        is_decreasing = f(theta, M) - f(theta, M-1) > log(L)

    if ((is_decreasing and f(theta, M) - log(- expm1(f(theta, M) - f(theta, M-1))) >= log(mpf(2)) + leps)
        or (not is_decreasing and f(theta, M) - log(1 - log(L)) >= log(mpf(2)) + leps)):
        raise ValueError("It is not possible to reach the stopping criterion with the given M.")

    log_terms[k] = f(theta, k)
    log_terms[k+1] = f(theta, k+1)
    k+=1

    while (log_terms[k] >= log_terms[k-1] or
           (is_decreasing and log_terms[k] - log(- expm1(log_terms[k] - log_terms[k-1])) >= log(mpf(2)) + leps) or 
           (not is_decreasing and log_terms[k] - log(1 - log(L)) >= log(mpf(2)) + leps) and
           k < M+initial_k):
        k+=1
        log_terms[k] = f(theta, k)
    
    log_sum = logsumexp(log_terms[initial_k:(k+1)])

    Bound1 = log_terms[k] - log(1 - log(L))
    Bound2 = log_terms[k] - log(- expm1(log_terms[k] - log_terms[k-1]))
    return (k-initial_k, logsumexp([log_sum, Bound1 - log(mpf(2)), Bound2 - log(mpf(2))]))
