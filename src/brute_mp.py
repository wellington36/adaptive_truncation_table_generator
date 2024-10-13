from utils.utils import logsumexp
from mpmath import mpf, log


def brute_mp(f, theta, M, initial_k):
    k = initial_k
    log_terms = [log(mpf(0))] * (M+initial_k+1)

    log_terms[k] = f(theta, k)
    log_terms[k+1] = f(theta, k+1)
    k+=1

    while (k < M+initial_k):
        k+=1
        log_terms[k] = f(theta, k)
    
    log_sum = logsumexp(log_terms)

    return (k-initial_k, log_sum)