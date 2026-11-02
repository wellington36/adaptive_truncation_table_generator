suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(purrr); library(tibble)
})

# Stable log PMF of total progeny J for NB offspring.
log_rj <- function(j, R0, k) {
  stopifnot(all(j >= 1), R0 > 0, k > 0)
  lgamma((k + 1) * j - 1) - lgamma(k * j) - lgamma(j + 1) +
    (j - 1) * log(R0 / k) - ((k + 1) * j - 1) * log1p(R0 / k)
}

# Adaptive enclosure. n indexes the last included term; evaluate next ratio
# from the exact pmf ratio. Returns an error if the cap is reached.
adaptive_rho <- function(R0, k, p_sent, eps_ll = 1e-8, m = 1L, max_n = 1e6L) {
  stopifnot(R0 > 0, k > 0, R0 < 1, p_sent > 0, p_sent < 1,
            eps_ll > 0, m >= 1, max_n >= 2)
  q <- 1 - p_sent
  n <- 1L
  log_r_n <- log_rj(1, R0, k)
  log_r_next <- log_rj(2, R0, k)

  log_t <- log_r_n + log(q)
  t_n <- exp(log_t)
  S_n <- t_n

  L <- q * R0 * ((k + 1) / (k + R0))^(k + 1)

  repeat {
    a_n <- exp(log(q) + log_r_next - log_r_n)

    if (!(a_n >= 0 && a_n < 1 && L >= a_n && L < 1))
      stop("Invalid ratio/bound; check arithmetic and model parameters.")

    lo <- S_n + t_n * a_n / (1 - a_n)
    hi <- S_n + t_n * L / (1 - L)

    if (!(lo >= 0 && lo <= hi && is.finite(lo) && !is.na(hi)))
      stop("Invalid rho interval ordering; check arithmetic and model parameters.")

    width_ll <- if (hi < 1)
      m * (log1p(-lo) - log1p(-hi)) else Inf

    if (width_ll <= eps_ll)
      return(list(lower = lo, upper = hi,
                  midpoint = (lo + hi)/2,
                  n = n, width_loglik = width_ll,
                  converged = TRUE))

    if (n >= max_n)
      stop("Maximum n reached before requested tolerance.")

    # Add the previously evaluated next term.
    tn1 <- exp(log_r_next + (n + 1) * log(q))

    S_n <- S_n + tn1
    t_n <- tn1
    n <- n + 1L

    # Advance the cached log-PMF values.
    log_r_n <- log_r_next
    log_r_next <- log_rj(n + 1, R0, k)
  }
}

# Fixed-cutoff approximation and observed log likelihood.
rho_fixed <- function(R0, k, p_sent, J) {
  j <- seq_len(J); q <- 1 - p_sent
  sum(exp(log_rj(j, R0, k) + j * log(q)))
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
  ll <- sum(log_rj(y, R0, k) + log1p(-q^y)) - m * log1p(-rho)
  list(logLik = ll, terms = terms, width = width)
}

# Detection-conditioned sampling. Build a sufficiently long table and check
# residual tail; increase J until the bound on omitted detected mass is small.
sample_detected <- function(m, R0, k, p_sent, J = 200000L, seed = 1L) {
  set.seed(seed); j <- seq_len(J); q <- 1 - p_sent
  logw <- log_rj(j, R0, k) + log1p(-q^j)
  w <- exp(logw - max(logw)); cdf <- cumsum(w) / sum(w)
  # Require stable truncation under a larger J in final production use.
  sample(j, size = m, replace = TRUE, prob = w)
}

# MLE helper: transform parameters to unconstrained scale; R0 logistic in (0,1),
# k and p_sent logistic transforms on practical bounded intervals to avoid
# optimizer excursions. Report/justify bounds before the final analysis.
fit_model <- function(y, method, J = 1000L, eps_ll = 1e-7, max_n = 1e6L,
                      starts = list(c(.5, 1, .2)), k_bounds = c(.05, 50),
                      p_bounds = c(.0005, .9995)) {
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
  #grid <- expand_grid(R0 = c(.3, .6, .9), k = c(.2, 1, 10),
  #                    p_sent = c(.01, .1, .5))
  grid <- expand_grid(R0 = c(.9), k = c(1.1), p_sent = c(.01, .1, .5))
  map_dfr(seq_len(nrow(grid)), function(g) {
    pars <- grid[g, ]
    map_dfr(seq_len(B), function(b) {
      y <- sample_detected(m, pars$R0, pars$k, pars$p_sent,
                           seed = seed + 100000L*g + b)
      # High-accuracy numerical reference: verify separately by tighter tolerances.
      ref <- fit_model(y, "adaptive", eps_ll = eps_ll/1000, max_n = 2e6L)
      methods <- c("adaptive", paste0("fixed_", fixed_J))
      map_dfr(methods, function(method) {
        tic <- proc.time()[["elapsed"]]
        if (method == "adaptive") fit <- fit_model(y, "adaptive", eps_ll = eps_ll)
        else fit <- fit_model(y, "fixed", J = as.integer(sub("fixed_", "", method)))
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
  aes(x = method, y = p95_abs_ll_error, color = method,
      group = 1)
) +
  geom_line(linewidth = .5) +
  geom_point(size = 2) +
  facet_grid(k ~ R0, labeller = label_both) +
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
  ggsave(file, p, width = 10, height = 7, dpi = 320); p
}

# ---------------- Experiment 2: published-data case study ----------------
run_published_data <- function(csv_path, size_column = "chain_size",
                               fixed_J = c(50L, 200L, 1000L), eps_ll = 1e-8) {
  dat <- read.csv(csv_path)
  if (!size_column %in% names(dat)) stop("CSV must contain column: ", size_column)
  y <- dat[[size_column]]
  if (!is.numeric(y) || anyNA(y) || any(y < 1) || any(y != floor(y)))
    stop("Chain sizes must be positive integer observations.")
  # Reproducibility record: report source, data version, preprocessing, and N.
  fits <- list(adaptive = fit_model(y, "adaptive", eps_ll = eps_ll,
                                     max_n = 2e6L))
  for (J in fixed_J) fits[[paste0("fixed_", J)]] <- fit_model(y, "fixed", J = J)
  bind_rows(lapply(names(fits), function(nm) {
    f <- fits[[nm]]; tibble(method = nm, R0 = as.numeric(f$R0), k = as.numeric(f$k),
       p_sent = as.numeric(f$p), logLik = as.numeric(f$logLik),
       convergence = as.integer(f$convergence))
  })) %>% mutate(delta_logLik = logLik - logLik[method == "adaptive"][1])
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

# ---------------- Results placeholders ----------------
# After running and validating the code, save output tables/figures and fill in:
# Experiment 1: [B, m, validated reference tolerance/cap], [key accuracy result],
# [runtime/term-count result], [convergence/coverage result].
# Experiment 2: [source/data preprocessing/N], [parameter estimates],
# [likelihood differences], [runtime/terms], [sensitivity conclusion].
