library(Rmpfr)

project_root_local <- if (exists("project_root", inherits = TRUE)) project_root else getwd()
source(file.path(project_root_local, "utils", "utils.R"))

sequential <- function(f, theta, M, eps, initial_k = 0L, prec = 64L) {
  prec <- as.integer(prec)
  M <- as.integer(M)
  initial_k <- as.integer(initial_k)
  eps <- mpfr(eps, prec)
  leps <- log(eps)

  term <- function(k) mpfr(f(theta, k), prec)

  # --- Compute the brute-force approximation with M terms ---
  ks <- seq.int(initial_k, initial_k + M)
  log_terms_brute <- do.call(c, lapply(ks, term))
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