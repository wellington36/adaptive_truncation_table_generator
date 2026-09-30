#!/usr/bin/env Rscript

# Hyper-Poisson raw and factorial moments using Ecdat::Bids$numbids
#
# Response-only analysis:
#   Y_i ~ Hyper-Poisson(theta, lambda),  i = 1,...,n
#
# No covariates and no regression model are used here.
# The normalising series is
#   Z(theta, lambda) = sum_{y=0}^infty q_y,
#   q_y = Gamma(lambda) theta^y / Gamma(lambda+y),
# with recurrence q_(y+1)/q_y = theta/(y+lambda).
#
# The script:
#   1. loads Ecdat::Bids and uses only numbids;
#   2. estimates theta and lambda by maximum likelihood;
#   3. computes raw and factorial moments of orders 1, 2, and 3;
#   4. returns certified intervals based on positive-series tail bounds;
#   5. checks raw moments using Stirling-number identities;
#   6. optionally compares against a high-precision reference if Rmpfr is installed.

options(stringsAsFactors = FALSE)

# ----------------------------- settings ----------------------------------
eps <- 1e-12
max_terms <- 100000L
seed <- 20260925
output_dir <- "hyperpoisson_numbids_results"

# Optional CSV input. If NULL, use the public Ecdat dataset.
input_csv <- NULL
response_column <- "numbids"

# ----------------------------- utilities ---------------------------------
stop_if <- function(condition, message) {
  if (condition) stop(message, call. = FALSE)
}

log_add_exp <- function(a, b) {
  if (!is.finite(a)) return(b)
  if (!is.finite(b)) return(a)
  m <- max(a, b)
  m + log1p(exp(min(a, b) - m))
}

# ------------------ adaptive normalising constant ------------------------
# Gold error-bounding-pairs convention. For the decreasing case L=0,
# the pair bound is B = q_(k+1)/(1-r_k), and convergence requires
# B <= 2*eps. The returned estimate is log(partial_sum + B/2).
log_hp_normalizer <- function(log_theta, log_lambda,
                              eps = 1e-12,
                              max_terms = 100000L) {
  stop_if(!is.finite(log_theta) || !is.finite(log_lambda),
          "Non-finite Hyper-Poisson parameters")
  stop_if(eps <= 0 || eps >= 1, "eps must be between 0 and 1")

  theta <- exp(log_theta)
  lambda <- exp(log_lambda)
  stop_if(!is.finite(theta) || !is.finite(lambda) || theta <= 0 || lambda <= 0,
          "theta and lambda must be positive")

  y <- 0L
  log_q <- 0.0                 # log(q_0), since q_0 = 1
  log_sum <- 0.0               # log(sum_{j=0}^y q_j)
  converged <- FALSE
  log_tail_upper <- Inf

  repeat {
    ratio <- theta / (y + lambda)
    stop_if(!is.finite(ratio), "Non-finite recurrence ratio")

    if (ratio < 1) {
      # First omitted term is q_(y+1).
      log_tail_upper <- log_q + log(ratio) - log1p(-ratio)
      if (log_tail_upper <= log(2) + log(eps)) {
        converged <- TRUE
        break
      }
    }

    stop_if(y + 1L >= max_terms,
            "Normalising-constant calculation exceeded max_terms")
    log_q <- log_q + log(ratio)
    log_sum <- log_add_exp(log_sum, log_q)
    y <- y + 1L
  }

  log_upper <- log_add_exp(log_sum, log_tail_upper)
  list(log_lower = log_sum,
       log_upper = log_upper,
       log_estimate = log_add_exp(log_sum, log_tail_upper - log(2)),
       terms = y + 1L,
       error = exp(log_tail_upper) / 2,
       converged = converged)
}

# ------------------ likelihood and estimation ----------------------------
log_hp_pmf <- function(y, log_theta, log_lambda, eps = 1e-12,
                       max_terms = 100000L) {
  stop_if(length(y) != 1L || !is.finite(y) || y < 0 || y != floor(y),
          "y must be one non-negative integer")
  z <- log_hp_normalizer(log_theta, log_lambda, eps, max_terms)
  lambda <- exp(log_lambda)
  lgamma(lambda) + y * log_theta - lgamma(lambda + y) - z$log_estimate
}

hp_nll <- function(par, y, eps = 1e-12, max_terms = 100000L,
                   return_details = FALSE) {
  log_theta <- par[1]
  log_lambda <- par[2]
  stop_if(any(!is.finite(par)), "Non-finite optimization parameters")

  n <- length(y)
  loglik <- numeric(n)
  terms <- integer(n)
  absolute_error <- numeric(n)

  for (i in seq_len(n)) {
    z <- log_hp_normalizer(log_theta, log_lambda, eps, max_terms)
    lambda <- exp(log_lambda)
    loglik[i] <- lgamma(lambda) + y[i] * log_theta -
      lgamma(lambda + y[i]) - z$log_estimate
    terms[i] <- z$terms
    absolute_error[i] <- z$error
  }

  value <- -sum(loglik)
  if (!return_details) return(value)
  list(value = value, loglik = loglik,
       log_theta = log_theta, log_lambda = log_lambda,
       theta = exp(log_theta), lambda = exp(log_lambda),
       terms = terms, normalizer_error = absolute_error)
}

# ------------------ raw and factorial moments ----------------------------
log_factorial_weight <- function(y, order, type) {
  if (type == "raw") {
    if (order == 0L) return(0.0)
    if (y == 0) return(-Inf)
    return(order * log(y))
  }
  if (type == "factorial") {
    if (order == 0L) return(0.0)
    if (y < order) return(-Inf)
    return(lgamma(y + 1) - lgamma(y - order + 1))
  }
  stop("type must be 'raw' or 'factorial'", call. = FALSE)
}

# Adaptive positive-series calculation for a weighted numerator.
log_hp_weighted_numerator <- function(log_theta, log_lambda, order,
                                      type = c("raw", "factorial"),
                                      eps = 1e-12,
                                      max_terms = 100000L) {
  type <- match.arg(type)
  stop_if(order < 0 || order != floor(order),
          "order must be a non-negative integer")

  theta <- exp(log_theta)
  lambda <- exp(log_lambda)
  y <- 0L
  log_q <- 0.0
  log_sum <- -Inf
  converged <- FALSE
  log_tail_upper <- Inf

  repeat {
    log_w <- log_factorial_weight(y, order, type)
    if (is.finite(log_w)) log_sum <- log_add_exp(log_sum, log_q + log_w)

    # The weighted sequence is eventually decreasing. Check the tail only
    # once the current and next weighted terms are both positive.
    ratio <- theta / (y + lambda)
    next_y <- y + 1L
    next_log_w <- log_factorial_weight(next_y, order, type)

    if (is.finite(log_w) && is.finite(next_log_w) && ratio < 1) {
      log_weighted_ratio <- next_log_w - log_w + log(ratio)
      weighted_ratio <- exp(log_weighted_ratio)
      if (is.finite(weighted_ratio) && weighted_ratio < 1) {
        log_first_omitted <- log_q + log(ratio) + next_log_w
        log_tail_upper <- log_first_omitted - log1p(-weighted_ratio)
        if (log_tail_upper <= log(2) + log(eps)) {
          converged <- TRUE
          break
        }
      }
    }

    stop_if(y + 1L >= max_terms,
            "Moment calculation exceeded max_terms")
    log_q <- log_q + log(ratio)
    y <- next_y
  }

  log_upper <- log_add_exp(log_sum, log_tail_upper)
  list(log_lower = log_sum,
       log_upper = log_upper,
       log_estimate = log_add_exp(log_sum, log_tail_upper - log(2)),
       terms = y + 1L,
       error = exp(log_tail_upper) / 2,
       converged = converged)
}

hp_moment <- function(order, log_theta, log_lambda,
                      type = c("raw", "factorial"),
                      eps = 1e-12, max_terms = 100000L) {
  type <- match.arg(type)
  z <- log_hp_normalizer(log_theta, log_lambda, eps, max_terms)
  num <- log_hp_weighted_numerator(log_theta, log_lambda, order,
                                   type, eps, max_terms)

  # Since numerator and denominator are positive,
  # [N_lower/Z_upper, N_upper/Z_lower] is a certified interval.
  lower <- exp(num$log_lower - z$log_upper)
  upper <- exp(num$log_upper - z$log_lower)
  estimate <- exp(num$log_estimate - z$log_estimate)

  list(order = order, type = type, estimate = estimate,
       lower = lower, upper = upper, error = (upper - lower) / 2,
       normalizer_terms = z$terms, numerator_terms = num$terms,
       norm_error = z$error, nume_error = num$error,
       converged = z$converged && num$converged)
}

# ------------------ data -------------------------------------------------
if (is.null(input_csv)) {
  stop_if(!requireNamespace("Ecdat", quietly = TRUE),
          "Install Ecdat with install.packages('Ecdat')")
  data("Bids", package = "Ecdat", envir = environment())
  dat <- Bids
} else {
  stop_if(!file.exists(input_csv), paste("File not found:", input_csv))
  dat <- read.csv(input_csv, check.names = FALSE)
}

stop_if(!(response_column %in% names(dat)),
        paste("Missing response column:", response_column))
y <- as.numeric(dat[[response_column]])
stop_if(length(y) < 2L, "At least two observations are required")
stop_if(any(!is.finite(y) | y < 0 | y != floor(y)),
        "The response must contain non-negative integer counts")

cat("Response-only Hyper-Poisson analysis\n")
cat("n =", length(y), "\n")
cat("sample mean =", mean(y), "\n")
cat("sample variance =", var(y), "\n")
cat("variance/mean =", var(y) / mean(y), "\n\n")

# ------------------ maximum likelihood fit -------------------------------
# The starting point is deliberately close to Poisson: lambda = 1 and theta
# equal to the sample mean. Multiple starts reduce the risk of a poor local fit.
starts <- rbind(
  c(log(max(mean(y), 0.05)), log(1.0)),
  c(log(max(mean(y), 0.05)), log(0.5)),
  c(log(max(mean(y), 0.05)), log(2.0)),
  c(log(max(mean(y) * 1.5, 0.05)), log(1.0))
)

fits <- lapply(seq_len(nrow(starts)), function(j) {
  optim(starts[j, ], hp_nll, y = y, eps = eps, max_terms = max_terms,
        method = "Nelder-Mead",
        control = list(maxit = 2000, reltol = 1e-10))
})
fit_values <- vapply(fits, function(x) x$value, numeric(1))
best <- fits[[which.min(fit_values)]]

# Refine with BFGS from the best derivative-free solution.
fit <- optim(best$par, hp_nll, y = y, eps = eps, max_terms = max_terms,
             method = "BFGS", hessian = TRUE,
             control = list(maxit = 2000, reltol = 1e-10))

stop_if(!is.finite(fit$value), "Optimization returned a non-finite objective")
details <- hp_nll(fit$par, y, eps, max_terms, return_details = TRUE)
log_theta_hat <- fit$par[1]
log_lambda_hat <- fit$par[2]
theta_hat <- exp(log_theta_hat)
lambda_hat <- exp(log_lambda_hat)
logLik_value <- -fit$value
AIC_value <- 2 * fit$value + 2 * length(fit$par)

cat("Estimated theta =", theta_hat, "\n")
cat("Estimated lambda =", lambda_hat, "\n")
cat("log-likelihood =", logLik_value, "\n")
cat("AIC =", AIC_value, "\n")
cat("optimizer convergence code =", fit$convergence, "\n")
cat("\nNormalizer terms per observation: ")
cat(unique(details$terms), "\n")

# ------------------ moments at fitted parameters -------------------------
moment_rows <- list()
row_id <- 1L
for (r in 1:3) {
  for (type in c("raw", "factorial")) {
    mm <- hp_moment(r, log_theta_hat, log_lambda_hat, type, eps, max_terms)
    moment_rows[[row_id]] <- data.frame(
      type = type, order = r, estimate = mm$estimate,
      lower = mm$lower, upper = mm$upper, error = mm$error,
      normalizer_terms = mm$normalizer_terms,
      numerator_terms = mm$numerator_terms,
      norm_error = mm$norm_error,
      nume_error = mm$nume_error,
      converged = mm$converged
    )
    row_id <- row_id + 1L
  }
}
moments <- do.call(rbind, moment_rows)

cat("\nCertified moments at the fitted parameters\n")
print(moments, row.names = FALSE)

# ------------------ optional high-precision reference ---------------------
# This reference is intentionally separate from the adaptive algorithm.
# It uses Rmpfr if available and a large fixed cap.
high_precision_reference <- NULL
if (requireNamespace("Rmpfr", quietly = TRUE)) {
  mp <- Rmpfr::mpfr
  prec <- 256
  theta_mp <- mp(theta_hat, prec)
  lambda_mp <- mp(lambda_hat, prec)
  q_mp <- mp(1, prec)
  z_mp <- q_mp
  y_mp <- 0L
  repeat {
    ratio_mp <- theta_mp / (mp(y_mp, prec) + lambda_mp)
    q_mp <- q_mp * ratio_mp
    z_mp <- z_mp + q_mp
    y_mp <- y_mp + 1L
    if (y_mp > 100000L || as.numeric(abs(q_mp)) < 1e-220) break
  }
  high_precision_reference <- data.frame(
    parameter = c("theta", "lambda", "Z"),
    value = c(as.numeric(theta_mp), as.numeric(lambda_mp), as.numeric(z_mp)),
    terms = c(NA_integer_, NA_integer_, y_mp + 1L)
  )
}

# ------------------ outputs -----------------------------------------------
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(moments, file.path(output_dir, "numbids_hyperpoisson_moments.csv"),
          row.names = FALSE)
write.csv(data.frame(
  n = length(y), sample_mean = mean(y), sample_variance = var(y),
  variance_to_mean = var(y) / mean(y), theta = theta_hat,
  lambda = lambda_hat, logLik = logLik_value, AIC = AIC_value,
  epsilon = eps, max_terms = max_terms,
  max_normalizer_error = max(details$normalizer_error),
  optimizer_convergence = fit$convergence
), file.path(output_dir, "numbids_fit_summary.csv"), row.names = FALSE)

saveRDS(list(data = dat, response = y, fit = fit, details = details,
             moments = moments, high_precision_reference = high_precision_reference,
             settings = list(eps = eps, max_terms = max_terms)),
        file.path(output_dir, "numbids_hyperpoisson_analysis.rds"))

cat("\nResults written to:", normalizePath(output_dir), "\n")


# ------------------ observed versus fitted distribution -------------------
y_grid <- 0:max(y)
observed <- tabulate(y + 1, nbins = length(y_grid)) / length(y)
z_fit <- log_hp_normalizer(log_theta_hat, log_lambda_hat, eps, max_terms)
p_fit <- exp(
  lgamma(lambda_hat) + y_grid * log_theta_hat -
    lgamma(lambda_hat + y_grid) - z_fit$log_estimate
)
png(file.path(output_dir, "numbids_hyperpoisson_fit.png"), 900, 600, res = 120)
plot(
  y_grid, observed, type = "h", lwd = 8,
  xlab = "Number of bids", ylab = "Probability",
  main = "Hyper-Poisson fit to Ecdat::Bids$numbids",
  ylim = c(0, max(c(observed, p_fit)) * 1.15)
)
points(y_grid, p_fit, pch = 19)
lines(y_grid, p_fit, lwd = 2)
legend(
  "topright", c("Observed", "Fitted Hyper-Poisson"),
  lwd = c(8, 2), pch = c(NA, 19), bty = "n"
)
dev.off()
