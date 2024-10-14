from utils.utils import logsumexp
from mpmath import mpf, log


def fixed_mp(f, theta, M, initial_k):
    '''Compute a fixed interations of the log-sum.
    
    Parameters:
    f                 (function): Log of the terms function
    theta (first parameter of f): Parameter of the function
    M                      (int): Maximum number of iterations
    initial_k              (int): Start of the sum

    Return: (iterations, approximation in log-scale)
    '''
    k = initial_k
    log_terms = [log(mpf(0))] * (M+initial_k+1)

    log_terms[k] = mpf(f(theta, k))
    log_terms[k+1] = mpf(f(theta, k+1))
    k+=1

    while (k < M+initial_k):
        k+=1
        log_terms[k] = mpf(f(theta, k))
    
    log_sum = logsumexp(log_terms[initial_k:(k+1)])

    return (k-initial_k, log_sum)