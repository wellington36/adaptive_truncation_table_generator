from utils.utils import logsumexp
from mpmath import mpf, log


def sum_to_threshold_mp(f, theta, M, L, eps, initial_k):
    ''' Applies the sum-to-threshold approach to obtain the approximate sum.
    
    Parameters:
    f                 (function): Log of the terms function
    theta (first parameter of f): Parameter of the function
    M                      (int): Maximum number of iterations
    L                    (float): Parameter of the test ratio
    eps                  (float): Error tolerance
    initial_k              (int): Start of the sum

    Return: (iterations, approximation in log-scale)
    '''
    k = initial_k
    L = mpf(L)
    eps = mpf(eps)
    M_bound = (L+1)/2
    leps = log(eps)
    log_terms = [log(mpf(0))] * (M+initial_k)

    def term(k):
        return mpf(f(theta, k))

    if (term(M) > leps + log(1-M_bound) - log(M_bound)):
        raise ValueError("It is not possible to reach the stopping criterion with the given M.")

    log_terms[k] = term(k)

    while (log_terms[k] > leps + log(1-M_bound) - log(M_bound) and k < M-1+initial_k):
        k+=1
        log_terms[k] = term(k)
    
    log_sum = logsumexp(log_terms[initial_k:(k+1)])

    return (k-initial_k, log_sum)