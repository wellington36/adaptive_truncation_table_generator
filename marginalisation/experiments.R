suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(progress)
})

# Stable log(exp(x) + exp(y)); handles -Inf inputs.
logspace_add <- function(x, y) {
  m <- pmax(x, y)
  out <- m + log1p(exp(-abs(x - y)))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}

# Stable log(sum(exp(x_i))).
log_sum_exp <- function(x) {
  if (!length(x)) return(-Inf)
  m <- max(x)
  if (is.infinite(m) && m < 0) return(-Inf)
  m + log(sum(exp(x - m)))
}

# Stable log PMF of total progeny J for NB offspring.
log_rj <- function(j, R0, k) {
  stopifnot(all(j >= 1), R0 > 0, k > 0)
  lgamma((k + 1) * j - 1) - lgamma(k * j) - lgamma(j + 1) +
    (j - 1) * log(R0 / k) - ((k + 1) * j - 1) * log1p(R0 / k)
}

# Adaptive enclosure. Consecutive terms are evaluated in vectorized batches to
# reduce R-loop overhead during optimizer calls; returned n is the first tested
# index meeting the tolerance. Tail and partial sums use log-domain arithmetic.
adaptive_rho <- function(R0, k, p_sent, eps_ll = 1e-8, m = 1L,
                         max_n = 1e6L, batch_size = 10L) {
  stopifnot(R0 > 0, k > 0, R0 < 1, p_sent > 0, p_sent < 1,
            eps_ll > 0, m >= 1, max_n >= 2, batch_size >= 1)
  q <- 1 - p_sent
  log_q <- log(q)
  n <- 1L
  log_r_n <- log_rj(1, R0, k)
  log_r_next <- log_rj(2, R0, k)
  log_t_n <- log_r_n + log_q
  log_S_n <- log_t_n
  L <- q * R0 * ((k + 1) / (k + R0))^(k + 1)
  log_L <- log(L)

  repeat {
    # Check the current index before advancing, so n is the first tested index
    # at which the requested likelihood-scale tolerance is achieved.
    a_n <- exp(log_q + log_r_next - log_r_n)
    if (!(a_n >= 0 && a_n < 1 && L >= a_n && L < 1))
      stop("Invalid ratio/bound; check arithmetic and model parameters.")
    log_tail_lower <- log_t_n + log(a_n) - log1p(-a_n)
    log_tail_upper <- log_t_n + log_L - log1p(-L)
    lo <- exp(logspace_add(log_S_n, log_tail_lower))
    hi <- exp(logspace_add(log_S_n, log_tail_upper))
    if (!(lo >= 0 && lo <= hi && is.finite(lo) && !is.na(hi)))
      stop("Invalid rho interval ordering; check arithmetic and model parameters.")
    width_ll <- if (hi < 1) m * (log1p(-lo) - log1p(-hi)) else Inf
    if (width_ll <= eps_ll)
      return(list(lower = lo, upper = hi, midpoint = (lo + hi)/2,
                  n = n, width_loglik = width_ll, converged = TRUE))
    if (n >= max_n)
      stop("Maximum n reached before requested tolerance.")

    # Evaluate a block of future adjacent terms and ratios in vectorized lgamma calls.
    b <- min(as.integer(batch_size), as.integer(max_n - n))
    js <- n + seq_len(b)
    log_r_batch <- log_rj(js, R0, k)
    log_r_next_batch <- log_rj(js + 1L, R0, k)
    log_t_batch <- log_r_batch + js * log_q
    a_batch <- exp(log_q + log_r_next_batch - log_r_batch)
    if (any(!is.finite(a_batch)) || any(a_batch < 0 | a_batch >= 1) ||
        any(L < a_batch))
      stop("Invalid ratio/bound; check arithmetic and model parameters.")

    # Compute partial sums through each candidate index with scaled cumulative sums.
    mx <- max(log_t_batch)
    log_additions <- mx + log(cumsum(exp(log_t_batch - mx)))
    log_S_candidates <- logspace_add(log_S_n, log_additions)
    log_tail_lower <- log_t_batch + log(a_batch) - log1p(-a_batch)
    log_tail_upper <- log_t_batch + log_L - log1p(-L)
    log_lo <- logspace_add(log_S_candidates, log_tail_lower)
    log_hi <- logspace_add(log_S_candidates, log_tail_upper)
    lo_batch <- exp(log_lo)
    hi_batch <- exp(log_hi)
    if (any(!is.finite(lo_batch)) || any(is.na(hi_batch)) ||
        any(lo_batch < 0) || any(lo_batch > hi_batch))
      stop("Invalid rho interval ordering; check arithmetic and model parameters.")
    width_batch <- rep(Inf, length(hi_batch))
    valid_hi <- hi_batch < 1
    width_batch[valid_hi] <- m * (log1p(-lo_batch[valid_hi]) -
                                  log1p(-hi_batch[valid_hi]))
    hit <- which(width_batch <= eps_ll)
    if (length(hit)) {
      i <- hit[1]
      return(list(lower = lo_batch[i], upper = hi_batch[i],
                  midpoint = (lo_batch[i] + hi_batch[i])/2,
                  n = as.integer(js[i]), width_loglik = width_batch[i],
                  converged = TRUE))
    }

    # Commit the whole block; all its terms are included in the new partial sum.
    n <- as.integer(js[b])
    log_S_n <- log_S_candidates[b]
    log_t_n <- log_t_batch[b]
    log_r_n <- log_r_batch[b]
    log_r_next <- log_r_next_batch[b]
  }
}

# Fixed-cutoff approximation and observed log likelihood.
# Positive series sums use log-sum-exp; the likelihood sums log contributions directly.
rho_fixed <- function(R0, k, p_sent, J) {
  stopifnot(R0 > 0, k > 0, p_sent > 0, p_sent < 1, J >= 1)
  j <- seq_len(J)
  log_terms <- log_rj(j, R0, k) + j * log1p(-p_sent)
  exp(log_sum_exp(log_terms))
}

loglik_observed <- function(y, R0, k, p_sent, method = c("adaptive", "fixed"),
                            J = 1000L, eps_ll = 1e-8, max_n = 1e6L) {
  method <- match.arg(method); m <- length(y)
  if (method == "adaptive") {
    z <- adaptive_rho(R0, k, p_sent, eps_ll = eps_ll, m = m, max_n = max_n)
    # midpoint is used for optimization; retain interval width for reporting.
    rho <- z$midpoint; terms <- z$n; width <- z$width_loglik
  } else {
    rho <- rho_fixed(R0, k, p_sent, J); terms <- J; width <- NA_real_
  }
  if (!is.finite(rho) || rho >= 1) return(list(logLik = -Inf, terms = terms, width = width))
  q <- 1 - p_sent
  # These are already logarithms of probabilities, so ordinary summation is correct.
  ll <- sum(log_rj(y, R0, k) + log1p(-q^y)) - m * log1p(-rho)
  list(logLik = ll, terms = terms, width = width)
}

# Detection-conditioned sampling. Build a sufficiently long table and check
# residual tail; increase J until the bound on omitted detected mass is small.
sample_detected <- function(m, R0, k, p_sent, J = 200000L, seed = 1L) {
  set.seed(seed); j <- seq_len(J); q <- 1 - p_sent
  logw <- log_rj(j, R0, k) + log1p(-q^j)
  # Normalize log weights using log-sum-exp, avoiding overflow/underflow.
  log_norm <- log_sum_exp(logw)
  w <- exp(logw - log_norm)
  # Require stable truncation under a larger J in final production use.
  sample(j, size = m, replace = TRUE, prob = w)
}

# MLE helper: transform parameters to unconstrained scale; R0 logistic in (0,1),
# k and p_sent logistic transforms on practical bounded intervals to avoid
# optimizer excursions. Report/justify bounds before the final analysis.
fit_model <- function(y, method, J = 1000L, eps_ll = 1e-7, max_n = 1e6L,
                      starts = list(c(.5, 1, .2)), k_bounds = c(.05, 50),
                      p_bounds = c(.0005, .9995), progress_bar = NULL) {
  decode <- function(z) c(R0 = min(1 - 1e-10, max(1e-10, plogis(z[1]))),
    k = k_bounds[1] + diff(k_bounds) * plogis(z[2]),
    p = p_bounds[1] + diff(p_bounds) * plogis(z[3]))

  objective <- function(z) {
    pars <- decode(z)
    -loglik_observed(y, pars[1], pars[2], pars[3], method, J, eps_ll, max_n)$logLik
  }

  fits <- lapply(starts, function(s) {
    z0 <- c(qlogis(s[1]), qlogis((s[2]-k_bounds[1])/diff(k_bounds)),
            qlogis((s[3]-p_bounds[1])/diff(p_bounds)))

    if (!is.null(progress_bar)) progress_bar$tick()

    optim(z0, objective, method = "BFGS", hessian = FALSE,
          control = list(maxit = 2000, reltol = 1e-10))
  })

  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  pars <- decode(best$par)

  # Recover numerical information at the final estimates.
  final_ll <- loglik_observed(
    y, pars[1], pars[2], pars[3],
    method = method, J = J, eps_ll = eps_ll, max_n = max_n
  )

  c(as.list(pars),
    logLik = final_ll$logLik,
    terms = final_ll$terms,
    width = final_ll$width,
    convergence = best$convergence,
    message = if (best$convergence == 0) "converged" else "check optimizer")
}

# ---------------- Experiment 1: simulation ----------------
run_simulation <- function(B = 200L, m = 500L, seed = 42L,
                           eps_ll = 1e-7, fixed_J = c(50L, 200L, 1000L)) {

  grid <- expand_grid(R0 = c(0.5, 0.7, 0.9), k = c(0.1, 1, 10),
                      p_sent = c(.01, .1, .5))

  total_optim <- nrow(grid) * B * (2L + length(fixed_J))

  pb <- progress_bar$new(
    format = "  Progress [:bar] :percent | :current/:total fits | Elapsed: :elapsed | ETA: :eta",
    total = total_optim,
    clear = TRUE,
    width = 65
  )

  map_dfr(seq_len(nrow(grid)), function(g) {
    pars <- grid[g, ]

    map_dfr(seq_len(B), function(b) {
      y <- sample_detected(m, pars$R0, pars$k, pars$p_sent,
                           seed = seed + 100000L*g + b)
      ref <- fit_model(y, "adaptive", eps_ll = eps_ll/1000,
                       max_n = 2e6L, progress_bar = pb)

      methods <- c("adaptive", paste0("fixed_", fixed_J))

      map_dfr(methods, function(method) {
        tic <- proc.time()[["elapsed"]]

        if (method == "adaptive") {
          fit <- fit_model(y, "adaptive", eps_ll = eps_ll,
                           progress_bar = pb)
        } else {
          fit <- fit_model(y, "fixed",
                           J = as.integer(sub("fixed_", "", method)),
                           progress_bar = pb)
        }

        elapsed <- proc.time()[["elapsed"]] - tic

        tibble(
          R0 = pars$R0, k = pars$k, p_sent = pars$p_sent,
          replicate = b,
          method = method,
          est_R0 = as.numeric(fit$R0),
          est_k = as.numeric(fit$k),
          est_p = as.numeric(fit$p),
          logLik = as.numeric(fit$logLik),
          ref_logLik = as.numeric(ref$logLik),
          abs_ll_error = abs(logLik - ref_logLik),
          terms = as.integer(fit$terms),
          width = as.numeric(fit$width),
          elapsed_sec = elapsed,
          convergence = as.integer(fit$convergence)
        )
      })
    })
  })
}

summarise_simulation <- function(res) {
  res %>%
    group_by(R0, k, p_sent, method) %>%
    summarise(
      n = n(),
      mean_abs_ll_error = mean(abs_ll_error, na.rm = TRUE),
      p95_abs_ll_error = quantile(abs_ll_error, .95, na.rm = TRUE),
      mean_R0_bias = mean(est_R0 - first(R0), na.rm = TRUE),
      RMSE_R0 = sqrt(mean((est_R0 - first(R0))^2, na.rm = TRUE)),
      median_terms = median(terms, na.rm = TRUE),
      max_terms = max(terms, na.rm = TRUE),
      median_runtime = median(elapsed_sec, na.rm = TRUE),
      convergence_rate = mean(convergence == 0, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_simulation <- function(summary_tbl, file = "experiment1_accuracy_cost.png") {
  p <- ggplot(
    summary_tbl,
    aes(x = method, y = p95_abs_ll_error, color = method, group = 1)
  ) +
    geom_line(linewidth = .5) +
    geom_point(size = 2) +
    facet_grid(
      rows = vars(k),
      cols = vars(R0, p_sent),
      labeller = label_both
    ) +
    scale_y_log10() +
    labs(
      x = NULL,
      y = "95th percentile absolute log-likelihood error (log scale)",
      title = "Likelihood accuracy across simulation regimes",
      subtitle = "Reference fits must be independently checked for numerical stability"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1)
    )

  ggsave(file, p, width = 14, height = 7, dpi = 320)
  p
}

# ---------------- Experiment 2: published-data case study ----------------
run_published_data <- function(csv_path, size_column = "chain_size",
                               fixed_J = c(50L, 200L, 1000L),
                               eps_ll = 1e-8, n_rep = 3L) {

  dat <- read.csv(csv_path)

  if (!size_column %in% names(dat))
    stop("CSV must contain column: ", size_column)

  y <- dat[[size_column]]

  if (!is.numeric(y) || anyNA(y) || any(y < 1) || any(y != floor(y)))
    stop("Chain sizes must be positive integer observations.")

  # Run each fit n_rep times and compute mean elapsed time
  benchmark_fit <- function(...) {
    times <- numeric(n_rep)
    fits <- vector("list", n_rep)

    for (i in seq_len(n_rep)) {
      timing <- system.time(
        fits[[i]] <- fit_model(...)
      )
      times[i] <- timing[["elapsed"]]
    }

    list(
      fit = fits[[1]],
      time = mean(times)
    )
  }

  fits <- list(
    adaptive = benchmark_fit(
      y, "adaptive", eps_ll = eps_ll, max_n = 2e6L
    ),
    ref = benchmark_fit(
      y, "adaptive", eps_ll = eps_ll / 1000, max_n = 2e6L
    )
  )

  for (J in fixed_J) {
    fits[[paste0("fixed_", J)]] <- benchmark_fit(
      y, "fixed", J = J
    )
  }

  bind_rows(lapply(names(fits), function(nm) {
    f <- fits[[nm]]$fit

    tibble(
      method = nm,
      R0 = as.numeric(f$R0),
      k = as.numeric(f$k),
      p_sent = as.numeric(f$p),
      terms = as.numeric(f$terms),
      logLik = as.numeric(f$logLik),
      convergence = as.integer(f$convergence),
      time = fits[[nm]]$time
    )
  })) %>%
    mutate(
      delta_logLik = logLik - logLik[method == "ref"][1]
    )
}

plot_case_study <- function(case_tbl, file = "experiment2_cutoff_sensitivity.png") {
  p <- ggplot(case_tbl, aes(x = method, y = delta_logLik, fill = method)) +
    geom_col(width = .7) + geom_hline(yintercept = 0, color = "grey35") +
    labs(x = NULL, y = "Log-likelihood difference from adaptive fit",
         title = "Sensitivity of published-data fit to truncation rule",
         subtitle = "Adaptive fit is a numerical reference, not a known true parameter") +
    theme_minimal(base_size = 11) + theme(legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(file, p, width = 8, height = 5, dpi = 320); p
}
