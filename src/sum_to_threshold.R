# Standalone R translation of src/sum_to_threshold_mp.py.

project_root_local <- if (exists("project_root", inherits = TRUE)) project_root else getwd()
source(file.path(project_root_local, "utils", "utils.R"))

sum_to_threshold <- function(f, theta, M, L, eps, initial_k = 0L, prec = 64L) {
  if (!requireNamespace("Rmpfr", quietly = TRUE)) {
    stop("Package 'Rmpfr' is required.")
  }

  prec <- as.integer(prec)
  M <- as.integer(M)
  initial_k <- as.integer(initial_k)
  L <- as.numeric(L)
  M_bound <- (L + 1) / 2
  leps <- log(mpfr(eps, prec))
  threshold <- leps + log1p(-M_bound) - log(M_bound)
  term <- function(k) mpfr(f(theta, k), prec)

  if (term(M) > threshold) {
    stop("It is not possible to reach the stopping criterion with the given M.")
  }
  k <- initial_k
  terms <- term(k)
  while (terms[[length(terms)]] > threshold && k < M - 1L + initial_k) {
    k <- k + 1L
    terms <- c(terms, term(k))
  }
  list(iterations = k - initial_k, value = logsumexp(terms, prec = prec))
}
