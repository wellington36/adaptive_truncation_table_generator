from src.sum_to_threshold_mp import sum_to_threshold_mp
from src.bounding_pairs_mp import bounding_pairs_mp

from mpmath import mp, mpf, log
from tabulate import tabulate


# function of each term in log scale
def f(a: float, n: int):
    return - mpf(2)*log(mpf(n)+1) - (mpf(n)+1)*log(a)

if __name__ == "__main__":
    mp.dps = 100

    a = [mpf(2), mpf('1.1'), mpf('1.01'), mpf('1.001'), mpf('1.0001'), mpf('1.00001')]
    L = [1/x for x in a]
    M = [5*10**6 for _ in range(len(a))]
    initial_k = 0

    stt_terms = [None] * len(a)
    bp_terms = [None] * len(a)
    percentage = [None] * len(a)
    for i in range(len(a)):
        bp_terms[i] = bounding_pairs_mp(f, (a[i]), M[i], L[i], mpf(2)**mpf(-52), initial_k)[0]
        stt_terms[i] = sum_to_threshold_mp(f, (a[i]), M[i], L[i], mpf(2)**mpf(-52), initial_k)[0]
        percentage[i] = bp_terms[i]/stt_terms[i] - 1
    
    data = []
    for idx, (stt, bp, perc) in enumerate(zip(stt_terms, bp_terms, percentage)):
        data.append([f"a={a[idx]}", stt, bp, perc])
    
    headers = ["", "STT", "BP", "Percent"]

    print(tabulate(data, headers, tablefmt="fancy_grid"))