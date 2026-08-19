from mpmath import mpf, log, expm1
from utils.utils import logsumexp


def bounding_pairs_mp(f, theta, M, L, eps, initial_k, bucket_size=20):
    k = initial_k
    L = mpf(L)
    eps = mpf(eps)
    leps = log(eps)
    converged = False

    def term(k):
        return mpf(f(theta, k))

    # Monotonicity test
    if L == 0:
        is_decreasing = True
    else:
        is_decreasing = term(M) - term(M - 1) > log(L)

    # M sufficiency
    if (
        (
            is_decreasing
            and term(M) - log(-expm1(term(M) - term(M - 1)))
            >= log(2) + leps
        )
        or (
            not is_decreasing
            and term(M) + log(L) - log(1 - L)
            >= log(2) + leps
        )
    ):
        raise ValueError(
            "It is not possible to reach the stopping criterion with the given M."
        )

    # Initialize with two consecutive terms
    prev_term = term(k)
    current_term = term(k + 1)
    k += 1

    log_terms = [prev_term, current_term]

    log_sum_total = log(mpf(0))

    while not converged and k < M + initial_k:
        k += 1
        prev_term = current_term
        current_term = term(k)

        log_terms.append(current_term)

        # Test convergence whenever the bucket is full.
        if len(log_terms) >= bucket_size:
            if (
                current_term <= prev_term
                and (
                    (
                        is_decreasing
                        and current_term
                        - log(-expm1(current_term - prev_term))
                        <= log(2) + leps
                    )
                    or (
                        not is_decreasing
                        and current_term
                        + log(L)
                        - log(1 - L)
                        <= log(2) + leps
                    )
                )
            ):
                converged = True

            log_sum_total = logsumexp(
                [log_sum_total] + log_terms
            )

            # Start a fresh bucket
            log_terms = []

    Bound1 = current_term + log(L) - log(1 - L)
    Bound2 = current_term - log(-expm1(current_term - prev_term))

    result = logsumexp(
        [log_sum_total, Bound1 - log(2), Bound2 - log(2)]
    )

    return (k - initial_k, result)
