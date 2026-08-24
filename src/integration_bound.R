project_root_local <- if (exists("project_root", inherits = TRUE)) project_root else getwd()
source(file.path(project_root_local, "utils", "utils.R"))

integration_bound <- function(f, g, theta, M, eps, initial_k = 0L, prec = 64L) {
  # Assuming that the series passes the integration test,
  # we obtain an approximation with guaranteed error, using Rmpfr.
  if (!requireNamespace("Rmpfr", quietly = TRUE)) {
    stop("Package 'Rmpfr' is required.")
  }
  prec <- as.integer(prec)
  M <- as.integer(M)
  initial_k <- as.integer(initial_k)
  leps <- log(mpfr(eps, prec))

  fterm <- function(k) mpfr(f(theta, k), prec)
  gterm <- function(k) mpfr(g(theta, k), prec)

  if (gterm(M) > log(2) + leps) {
    stop("It is not possible to reach the stopping criterion with the given M.")
  }

  k <- initial_k
  terms <- mpfr(numeric(0), prec)

  while (gterm(k) > log(2) + leps) {
    terms <- c(terms, fterm(k))
    k <- k + 1L
    if (k > M) stop("integration bound exceeded M")
  }

  log_sum <- logsumexp(terms, prec = prec)
  result <- logsumexp(c(log_sum, gterm(k) - log(2), gterm(k - 1L) - log(2)), prec = prec)

  list(iterations = k - initial_k, value = result)
}