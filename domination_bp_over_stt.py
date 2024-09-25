from src.sum_to_threshold_mp import sum_to_threshold_mp
from src.bounding_pairs_mp import bounding_pairs_mp
from src.sequential_mp import sequential_mp
from src.brute_mp import brute_mp

# pip install /path/to/local/clone
import pybind_stan_fns as psf

from rpy2.robjects.packages import importr
from mpmath import mp, mpf, log, exp
from utils.utils import logdiffexp
from tabulate import tabulate
from math import lgamma
import math

def f(a: float, n: int):
    return - mpf(2)*log(mpf(n)+1) - (mpf(n)+1)*log(a)

if __name__ == "__main__":
    mp.dps = 200

    a = [mpf(100), mpf(10), mpf(2), mpf('1.1'), mpf('1.01'), mpf('1.001')]
    L = [1/x for x in a]
    M = [3*10**6 for _ in range(len(a))]

    stt_terms = [None] * len(a)
    bp_terms = [None] * len(a)
    percentage = [None] * len(a)
    for i in range(len(a)):
        bp_terms[i] = bounding_pairs_mp(f, (a[i]), M[i], L[i], mpf(2)**mpf(-52), initial_k=2)[0]
        stt_terms[i] = sum_to_threshold_mp(f, (a[i]), M[i], L[i], mpf(2)**mpf(-52), initial_k=2)[0]
        percentage[i] = bp_terms[i]/stt_terms[i] - 1
    
    data = []
    for idx, (stt, bp, perc) in enumerate(zip(stt_terms, bp_terms, percentage), start=1):
        data.append([f"a={a[idx-1]}", stt, bp, perc])
    
    headers = ["", "STT", "BP", "Percent"]

    print(tabulate(data, headers, tablefmt="fancy_grid"))