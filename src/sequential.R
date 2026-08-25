project_root_local <- if (exists("project_root", inherits = TRUE)) project_root else getwd()
source(file.path(project_root_local, "utils", "utils.R"))

sequential <- function(f, theta, M, eps, initial_k = 0L, prec = 64L) {
  # Approximate the sum with M terms and then evaluate iterations until
  # reaching the error with the approximated sum, using Rmpfr.
  if (!requireNamespace("Rmpfr", quietly = TRUE)) {
    stop("Package 'Rmpfr' is required.")
  }
  prec <- as.integer(prec)
  M <- as.integer(M)
  initial_k <- as.integer(initial_k)
  eps <- mpfr(eps, prec)
  leps <- log(eps)

  term <- function(k) mpfr(f(theta, k), prec)

  # --- Compute the brute-force approximation with M terms ---
  log_terms_brute <- mpfr(rep(-Inf, M + 1L), prec)

  for (k in seq.int(initial_k, M + initial_k)) {
    log_terms_brute[k - initial_k + 1L] <- term(k)
  }

  log_sum_brute <- logsumexp(log_terms_brute, prec = prec)

  # --- Compute the terms until we reach the approximation ---
  k <- initial_k
  log_sum <- logsumexp(log_terms_brute[1L:2L], prec = prec)
  k <- k + 1L

  repeat {
    if (log_sum >= log_sum_brute) break
    if (logdiffexp(log_sum_brute, log_sum, prec = prec) < leps) break
    if (k >= M + initial_k) break

    k <- k + 1L
    log_sum <- logsumexp(c(log_sum, log_terms_brute[k - initial_k + 1L]), prec = prec)
  }

  list(iterations = k - initial_k, value = log_sum)
}