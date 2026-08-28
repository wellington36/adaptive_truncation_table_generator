library(Rmpfr)

project_root_local <- if (exists("project_root", inherits = TRUE)) project_root else getwd()
source(file.path(project_root_local, "utils", "utils.R"))

bounding_pairs <- function(f, theta, M, L, eps, initial_k = 0L, bucket_size = 20L, prec = 64L) {
  prec <- as.integer(prec)
  k <- as.integer(initial_k)
  M <- as.integer(M)
  L <- as.numeric(L)
  leps <- log(mpfr(eps, prec))
  converged <- FALSE
  term <- function(j) mpfr(f(theta, j), prec)

  is_decreasing <- if (L == 0) TRUE else term(M) - term(M - 1L) > log(L)

  if ((is_decreasing && term(M) - log(-expm1(term(M) - term(M - 1L))) >= log(2) + leps) ||
      (!is_decreasing && term(M) + log(L) - log1p(-L) >= log(2) + leps)) {
    stop("It is not possible to reach the stopping criterion with the given M.")
  }

  prev_term <- term(k)
  current_term <- term(k + 1L)
  k <- k + 1L
  log_terms <- c(prev_term, current_term)
  log_sum_total <- mpfr(-Inf, prec)

  while (!converged && k < M + initial_k) {
    k <- k + 1L
    prev_term <- current_term
    current_term <- term(k)
    log_terms <- c(log_terms, current_term)
    if (length(log_terms) >= bucket_size) {
      if (current_term <= prev_term &&
          ((is_decreasing && current_term - log(-expm1(current_term - prev_term)) <= log(2) + leps) ||
           (!is_decreasing && current_term + log(L) - log1p(-L) <= log(2) + leps))) {
        converged <- TRUE
      }

      log_sum_total <- logsumexp(c(log_terms, log_sum_total), prec=prec)
      log_terms <- mpfr(numeric(0), prec)
    }
  }
  
  bound1 <- current_term + log(L) - log1p(-L)
  bound2 <- current_term - log(-expm1(current_term - prev_term))
  result <- logsumexp(c(log_sum_total, bound1 - log(2), bound2 - log(2)), prec=prec)
  list(iterations = k - initial_k, value = result)
}
