from utils.utils import logsumexp
from mpmath import mpf, log

def integration_bound_mp(f, g, theta, M, eps, initial_k):
    ''' Assuming that the series passes the integration test, we obtain an approximation with guaranteed error.
    
    Parameters:
    f                 (function): Log of the terms function
    g                 (function): Integral of n to infinity of f
    theta (first parameter of f): Parameter of the function
    M                      (int): Maximum number of iterations
    eps                  (float): Error tolerance
    initial_k              (int): Start of the sum

    Return: (iterations, approximation in log-scale)
    '''
    k = initial_k
    eps = mpf(eps)
    leps = log(eps)
    log_terms = [log(mpf(0))] * (M + initial_k)

    if (mpf(g(theta, M)) > log(mpf(2)) + leps):
        raise ValueError("It is not possible to reach the stopping criterion with the given M.")

    log_terms[k] = mpf(f(theta, k))
    k+=1

    while (mpf(g(theta, k)) > log(mpf(2)) + leps):
        log_terms[k] = mpf(f(theta, k))
        k+=1
    
    log_sum = logsumexp(log_terms[initial_k:(k+1)])

    return (k-initial_k, logsumexp([log_sum, g(theta, k) - log(mpf(2)), g(theta, k-1) - log(mpf(2))]))

