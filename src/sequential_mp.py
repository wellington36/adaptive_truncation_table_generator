from utils.utils import logsumexp, logdiffexp
from mpmath import mpf, log


def sequential_mp(f, theta, M, eps, initial_k):
    '''Approximate the sum with M terms and then evaluate iterations until reaching the error with the approximated sum.
    
    Parameters:
    f                 (function): Log of the terms function
    theta (first parameter of f): Parameter of the function
    M                      (int): Maximum number of iterations
    eps                  (float): Error tolerance
    initial_k              (int): Start of the sum

    Return: (iterations, approximation in log-scale)
    '''
    k = initial_k
    eps = mpf(eps)
    leps = log(eps)
    log_terms_brute = [log(mpf(0))] * (M+initial_k+1)

    # Aux function
    def term(k):
        return mpf(f(theta, k))

    # Compute the approximation
    log_terms_brute[k] = term(k)
    log_terms_brute[k+1] = term(k+1)
    k+=1

    while (k < M+initial_k):
        k+=1
        log_terms_brute[k] = term(k)
    
    log_sum_brute = logsumexp(log_terms_brute[initial_k:(k+1)])

    # Compute the terms until we reach the approximation
    k = initial_k
    log_terms = [log(mpf(0))] * (M+initial_k+1)

    log_terms[k] = term(k)
    log_terms[k+1] = term(k+1)
    k+=1
    
    log_sum = logsumexp(log_terms[initial_k:(k+1)])

    while (logdiffexp(log_sum_brute, log_sum) >= leps):
        k+=1
        log_terms[k] = term(k)
        log_sum = logsumexp([log_sum, log_terms[k]])
    
    return (k-initial_k, log_sum)