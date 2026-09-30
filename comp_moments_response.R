#!/usr/bin/env Rscript

# Response-only Conway-Maxwell-Poisson experiment
# Uses MASS::quine$Days and no regression covariates.
#
# p(y | lambda, nu) = lambda^y / (y!)^nu / Z(lambda, nu)
# q_y = exp(y log(lambda) - nu lgamma(y+1))
# q_(y+1)/q_y = lambda/(y+1)^nu

options(stringsAsFactors = FALSE)

if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Install MASS with install.packages('MASS')", call. = FALSE)
}
data("quine", package = "MASS")
y <- as.integer(quine$Days)
stopifnot(all(is.finite(y)), all(y >= 0L), all(y == floor(y)))

eps <- 1e-12
max_terms <- 100000L
output_dir <- "comp_moments_results"

log_add_exp <- function(a, b) {
  if (!is.finite(a)) return(b)
  if (!is.finite(b)) return(a)
  m <- max(a, b)
  m + log1p(exp(min(a, b) - m))
}

log_midpoint <- function(a, b) {
  a + log1p(exp(b - a)) - log(2)
}

# Certified lower/upper log bounds for Z(lambda, nu).
log_comp_normalizer <- function(log_lambda, log_nu,
                                eps = 1e-12,
                                max_terms = 1000000L) {
  lambda <- exp(log_lambda)
  nu <- exp(log_nu)
  stopifnot(lambda > 0, nu > 0, eps > 0, eps < 1)

  y <- 0L
  log_q <- 0
  log_sum <- 0
  converged <- FALSE
  log_tail_upper <- Inf

  repeat {
    ratio <- lambda / (y + 1)^nu
    if (ratio < 1) {
      log_tail_upper <- log_q + log(ratio) - log1p(-ratio)
      if (log_tail_upper <= log(2) + log(eps)) {
        converged <- TRUE
        break
      }
    }
    if (y + 1L >= max_terms) {
      stop("COMP normalizer exceeded max_terms", call. = FALSE)
    }
    log_q <- log_q + log(ratio)
    log_sum <- log_add_exp(log_sum, log_q)
    y <- y + 1L
  }

  log_upper <- log_add_exp(log_sum, log_tail_upper)
  list(log_lower = log_sum, log_upper = log_upper,
       log_estimate = log_add_exp(log_sum, log_tail_upper - log(2)),
       terms = y + 1L,
       error = exp(log_tail_upper)/2,
       converged = converged)
}

log_comp_lpmf <- function(y, log_lambda, log_nu,
                          eps = 1e-12, max_terms = 1000000L) {
  z <- log_comp_normalizer(log_lambda, log_nu, eps, max_terms)
  y * log_lambda - exp(log_nu) * lgamma(y + 1) - z$log_estimate
}

comp_nll <- function(par, y, eps = 1e-12,
                     max_terms = 1000000L, return_details = FALSE) {
  ll <- vapply(y, log_comp_lpmf, numeric(1),
               log_lambda = par[1], log_nu = par[2],
               eps = eps, max_terms = max_terms)
  if (!return_details) return(-sum(ll))
  list(value = -sum(ll), loglik = ll,
       lambda = exp(par[1]), nu = exp(par[2]))
}

log_weight <- function(y, order, type) {
  if (type == "raw") {
    if (order == 0L) return(0)
    if (y == 0L) return(-Inf)
    return(order * log(y))
  }
  if (type == "factorial") {
    if (order == 0L) return(0)
    if (y < order) return(-Inf)
    return(lgamma(y + 1) - lgamma(y - order + 1))
  }
  stop("type must be raw or factorial", call. = FALSE)
}

log_comp_weighted_numerator <- function(log_lambda, log_nu, order,
                                        type = c("raw", "factorial"),
                                        eps = 1e-12,
                                        max_terms = 1000000L) {
  type <- match.arg(type)
  lambda <- exp(log_lambda)
  nu <- exp(log_nu)
  y <- 0L
  log_q <- 0
  log_sum <- -Inf
  converged <- FALSE
  log_tail_upper <- Inf

  repeat {
    lw <- log_weight(y, order, type)
    if (is.finite(lw)) log_sum <- log_add_exp(log_sum, log_q + lw)

    next_y <- y + 1L
    ratio <- lambda / (y + 1)^nu
    next_lw <- log_weight(next_y, order, type)

    if (is.finite(lw) && is.finite(next_lw) && ratio < 1) {
      weighted_ratio <- exp(next_lw - lw) * ratio
      if (is.finite(weighted_ratio) && weighted_ratio < 1) {
        log_first_omitted <- log_q + log(ratio) + next_lw
        log_tail_upper <- log_first_omitted - log1p(-weighted_ratio)
        if (log_tail_upper <= log(2) + log(eps)) {
          converged <- TRUE
          break
        }
      }
    }

    if (y + 1L >= max_terms) {
      stop("COMP weighted numerator exceeded max_terms", call. = FALSE)
    }
    log_q <- log_q + log(ratio)
    y <- next_y
  }

  log_upper <- log_add_exp(log_sum, log_tail_upper)
  list(log_lower = log_sum, log_upper = log_upper,
       log_estimate = log_add_exp(log_sum, log_tail_upper - log(2)),
       terms = y + 1L,
       error = exp(log_tail_upper)/2,
       converged = converged)
}

comp_moment <- function(order, log_lambda, log_nu,
                        type = c("raw", "factorial"),
                        eps = 1e-12, max_terms = 1000000L) {
  type <- match.arg(type)
  z <- log_comp_normalizer(log_lambda, log_nu, eps, max_terms)
  n <- log_comp_weighted_numerator(log_lambda, log_nu, order, type,
                                   eps, max_terms)
  lower <- exp(n$log_lower - z$log_upper)
  upper <- exp(n$log_upper - z$log_lower)
  list(order = order, type = type,
       estimate = exp(n$log_estimate - z$log_estimate),
       lower = lower, upper = upper, error = (upper - lower)/2,
       normalizer_terms = z$terms, numerator_terms = n$terms,
       norm_error = z$error,
       nume_error = n$error,
       converged = z$converged && n$converged)
}

# Fit response-only COMP model by multiple-start optimization.
starts <- rbind(
  c(log(max(mean(y), 0.05)), log(1)),
  c(log(max(mean(y), 0.05)), log(0.5)),
  c(log(max(mean(y), 0.05)), log(2)),
  c(log(max(mean(y) * 1.5, 0.05)), log(1))
)
fits <- lapply(seq_len(nrow(starts)), function(i) {
  optim(starts[i, ], comp_nll, y = y, eps = eps,
        max_terms = max_terms, method = "Nelder-Mead",
        control = list(maxit = 2000, reltol = 1e-10))
})
best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
fit <- optim(best$par, comp_nll, y = y, eps = eps,
             max_terms = max_terms, method = "BFGS", hessian = TRUE,
             control = list(maxit = 2000, reltol = 1e-10))

lambda_hat <- exp(fit$par[1])
nu_hat <- exp(fit$par[2])
cat("Response-only COMP analysis\n")
cat("n =", length(y), "sample mean =", mean(y),
    "sample variance =", var(y), "\n")
cat("lambda =", lambda_hat, "nu =", nu_hat, "\n")
cat("log-likelihood =", -fit$value, "AIC =", 2 * fit$value + 4, "\n")

rows <- list()
k <- 1L
for (r in 1:3) {
  for (type in c("raw", "factorial")) {
    m <- comp_moment(r, fit$par[1], fit$par[2], type, eps, max_terms)
    rows[[k]] <- data.frame(type = type, order = r,
      estimate = m$estimate, lower = m$lower, upper = m$upper,
      error = m$error, normalizer_terms = m$normalizer_terms,
      numerator_terms = m$numerator_terms,
      norm_error = m$norm_error,
      nume_error = m$nume_error,
      converged = m$converged)
    k <- k + 1L
  }
}
moments <- do.call(rbind, rows)
print(moments, row.names = FALSE)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(moments, file.path(output_dir, "comp_days_moments.csv"), row.names = FALSE)
write.csv(data.frame(n = length(y), sample_mean = mean(y),
  sample_variance = var(y), lambda = lambda_hat, nu = nu_hat,
  logLik = -fit$value, AIC = 2 * fit$value + 4, epsilon = eps,
  max_terms = max_terms),
  file.path(output_dir, "comp_days_fit_summary.csv"), row.names = FALSE)
saveRDS(list(y = y, fit = fit, moments = moments,
             settings = list(eps = eps, max_terms = max_terms)),
        file.path(output_dir, "comp_days_analysis.rds"))
cat("Results written to", normalizePath(output_dir), "\n")


# Observed vs fitted COMP
y_grid <- 0:max(y)
observed <- tabulate(y + 1, nbins = length(y_grid)) / length(y)

z_fit <- log_comp_normalizer(
  fit$par[1], fit$par[2], eps, max_terms
)

p_fit <- exp(
  y_grid * fit$par[1] -
  exp(fit$par[2]) * lgamma(y_grid + 1) -
  z_fit$log_estimate
)

png(file.path(output_dir, "comp_days_fit.png"), 900, 600, res = 120)

plot(
  y_grid, observed, type = "h", lwd = 8,
  xlab = "Days", ylab = "Probability",
  main = "Conway-Maxwell-Poisson fit to quine$Days",
  ylim = c(0, max(c(observed, p_fit)) * 1.15)
)
points(y_grid, p_fit, pch = 19)
lines(y_grid, p_fit, lwd = 2)

legend(
  "topright", c("Observed", "Fitted COMP"),
  lwd = c(8, 2), pch = c(NA, 19), bty = "n"
)

dev.off()