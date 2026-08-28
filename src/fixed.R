library(Rmpfr)

project_root_local <- if (exists("project_root", inherits = TRUE)) project_root else getwd()
source(file.path(project_root_local, "utils", "utils.R"))

fixed <- function(f, theta, M, initial_k = 0L, prec = 64L) {
  # Compute a fixed number of iterations of the log-sum using Rmpfr.

  prec <- as.integer(prec)
  initial_k <- as.integer(initial_k)
  M <- as.integer(M)

  if (M < 1L || initial_k < 0L)
    stop("M must be positive and initial_k non-negative")

  term <- function(j) mpfr(f(theta, j), prec)

  ks <- seq.int(initial_k, initial_k + M)
  log_terms <- do.call(c, lapply(ks, term))
  log_sum <- logsumexp(log_terms, prec = prec)

  list(iterations = M, value = log_sum)
}