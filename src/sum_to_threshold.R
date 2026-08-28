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
  log_M_bound <- log(mpfr(M_bound, prec))
  leps <- log(mpfr(eps, prec))
  threshold <- leps + log1p(-M_bound) - log_M_bound
  term <- function(k) mpfr(f(theta, k), prec)

  if (term(M) > threshold) {
    stop("It is not possible to reach the stopping criterion with the given M.")
  }

  k <- initial_k
  terms <- term(k)
  prev_term <- terms[[length(terms)]]

  ratio_ok <- FALSE  # no previous pair yet, so nothing to violate

  while ((prev_term > threshold || !ratio_ok) && k < M - 1L + initial_k) {
    k <- k + 1L
    new_term <- term(k)

    ratio_ok <- (new_term - prev_term) <= log_M_bound

    terms <- c(terms, new_term)
    prev_term <- new_term
  }

  list(
    iterations = k - initial_k,
    value = logsumexp(terms, prec = prec)
  )
}