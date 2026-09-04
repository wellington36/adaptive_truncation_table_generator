#############################################################################
## Adaptive truncation for the ZICOM-Poisson normalising constant
## Reproducing / stress-testing Barriga & Louzada (2014, Stat. Methodology)
## "The zero-inflated Conway-Maxwell-Poisson distribution: Bayesian
##  inference, regression modeling and influence diagnostic"
## Section 6.2: Apple cultivar (Trajan) data.
##
## Goal: compare the paper's *fixed* truncation at K = 100 terms of
##   S(mu, nu) = sum_{s=0}^Inf (mu^s / s!)^nu
## against two adaptive, tolerance-guaranteed truncation schemes:
##   (1) Sum-to-threshold
##   (2) Error-bounding pairs
## inside a Metropolis-Hastings sampler for the ZICOM-Poisson regression
## model, using the same priors and MH design as the original paper
## (independent N(0, 10^2) priors, random-walk MH).
##
## Data: Ridout, Hinde & Demetrio (1998) "Models for count data with many
## zeros" -- 270 micropropagated shoots of the columnar apple cultivar
## Trajan.
#############################################################################

set.seed(2024)

## ---------------------------------------------------------------------
## 0. Data
## ---------------------------------------------------------------------
if (!requireNamespace("agridat", quietly = TRUE)) {
  stop("Please run install.packages('agridat') first (CRAN).")
}
data(ridout.appleshoots, package = "agridat")
dat <- ridout.appleshoots
dat$photo01 <- as.numeric(dat$photo == 16)  # 0 = 8h, 1 = 16h

cat("\nphoto / photo01 cross-tab (confirm this is 2 clean groups):\n")
print(table(raw = dat$photo, photo01 = dat$photo01))
cat(sprintf("n = %d, %% zero = %.1f%%\n", nrow(dat), 100 * mean(dat$roots == 0)))

## Raw empirical zero rates by photoperiod, for the independent sanity
## check later (from the paper's own Sec. 6.2 description / directly
## from the data -- NOT from Table 4, so this isn't circular).
raw_zero_8h  <- mean(dat$roots[dat$photo01 == 0] == 0)
raw_zero_16h <- mean(dat$roots[dat$photo01 == 1] == 0)
cat(sprintf("Raw P(roots=0) at 8h = %.4f, at 16h = %.4f\n", raw_zero_8h, raw_zero_16h))

## ---------------------------------------------------------------------
## 1. Truncation strategies
## ---------------------------------------------------------------------

## (a) Fixed truncation -- what the original paper used (K = 100, i.e.
##     summing the first 101 terms, "truncated at the 101st term").
S_fixed <- function(mu, nu, K = 100) {
  log_mu <- log(mu); logA <- 0; S <- 1; n_terms <- 1L
  for (n in 1:K) {
    logA <- logA + nu * (log_mu - log(n))
    a <- exp(logA)
    if (!is.finite(a)) return(list(S = NA_real_, n_terms = n_terms))
    S <- S + a
    n_terms <- n_terms + 1L
  }
  list(S = S, n_terms = n_terms)
}

## (b) Sum-to-threshold (Approach 1): stop once a_n < eps*(1-M)/M AND the
##     ratio has dropped below M. Here L = 0 (CMP-type ratio -> 0), so we
##     take M = 0.5 as a safe interior point of (L, 1).
S_threshold <- function(mu, nu, eps = 2.2e-10, M = 0.5, n_max = 2e5) {
  log_mu <- log(mu); logA_prev <- 0; S <- 1; n <- 0L
  thresh <- eps * (1 - M) / M
  repeat {
    n <- n + 1L
    log_r <- nu * (log_mu - log(n))
    logA_cur <- logA_prev + log_r
    a_cur <- exp(logA_cur)
    if (!is.finite(a_cur)) return(list(S = NA_real_, n_terms = n))
    S <- S + a_cur
    r <- exp(log_r)
    if (a_cur < thresh && r < M) break
    if (n >= n_max) break
    logA_prev <- logA_cur
  }
  list(S = S, n_terms = n + 1L)
}

## (c) Error-bounding pairs (Approach 2): stop once the ratio is
##     decreasing (r < 1) and the bound a_n/(1-r) drops below 2*eps;
##     return the midpoint-corrected sum, per Eq. (adaptSum) in the
##     truncation paper.
S_bounding <- function(mu, nu, eps = 2.2e-10, n_max = 2e5) {
  log_mu <- log(mu); logA_prev <- 0; S <- 1; n <- 0L; a_cur <- 1; r <- Inf
  repeat {
    n <- n + 1L
    log_r <- nu * (log_mu - log(n))
    logA_cur <- logA_prev + log_r
    a_cur <- exp(logA_cur)
    if (!is.finite(a_cur)) return(list(S = NA_real_, n_terms = n))
    S <- S + a_cur
    r <- exp(log_r)
    if (r < 1 && a_cur / (1 - r) < 2 * eps) break
    if (n >= n_max) break
    logA_prev <- logA_cur
  }
  correction <- if (r < 1) (a_cur / 2) * (1 / (1 - r)) else 0
  list(S = S + correction, n_terms = n + 1L)
}

## ---------------------------------------------------------------------
## 2. ZICOM-Poisson log-likelihood / log-prior / log-posterior
##    (Barriga & Louzada, Eqs. (4)-(7); Section 6.2 link functions)
## ---------------------------------------------------------------------
loglik <- function(theta, dat, Sfun, ...) {
  a0 <- theta[1]; a1 <- theta[2]; a2 <- theta[3]
  b0 <- theta[4]; b1 <- theta[5]; b2 <- theta[6]
  g0 <- theta[7]; g1 <- theta[8]; g2 <- theta[9]

  p   <- plogis(a0 + a1 * dat$bap + a2 * dat$photo01)
  mu  <- exp(b0 + b1 * dat$bap + b2 * dat$photo01)
  nu  <- exp(g0 + g1 * dat$bap + g2 * dat$photo01)

  ll <- 0; total_terms <- 0L
  for (i in seq_len(nrow(dat))) {
    res <- Sfun(mu[i], nu[i], ...)
    S <- res$S
    total_terms <- total_terms + res$n_terms
    if (!is.finite(S) || S < 1) return(list(ll = -Inf, n_terms = total_terms))
    if (dat$roots[i] == 0) {
      val <- p[i] + (1 - p[i]) / S
      if (!is.finite(val) || val <= 0) return(list(ll = -Inf, n_terms = total_terms))
      ll <- ll + log(val)
    } else {
      y <- dat$roots[i]
      ll <- ll + log(1 - p[i]) + nu[i] * (y * log(mu[i]) - lgamma(y + 1)) - log(S)
    }
  }
  if (!is.finite(ll)) return(list(ll = -Inf, n_terms = total_terms))
  list(ll = ll, n_terms = total_terms)
}

logprior <- function(theta, sd = 10) sum(dnorm(theta, 0, sd, log = TRUE))

logpost <- function(theta, dat, Sfun, ...) {
  lp <- logprior(theta)
  if (!is.finite(lp)) return(list(lp = -Inf, n_terms = 0L))
  out <- loglik(theta, dat, Sfun, ...)
  list(lp = lp + out$ll, n_terms = out$n_terms)
}

## ---------------------------------------------------------------------
## 3. Sampler: joint move, DIAGONAL-only proposal covariance, auto-tuned
##    via windowed Robbins-Monro adaptation (no hand-coded step sizes,
##    no off-diagonal correlation estimation).
## ---------------------------------------------------------------------
tune_diagonal_scale <- function(dat, Sfun, theta0, n_pilot = 3000,
                                 init_sd = rep(0.05, length(theta0)),
                                 adapt_every = 100, window = 500,
                                 target_accept = 0.234, ...) {
  d <- length(theta0)
  theta <- theta0; cur <- logpost(theta, dat, Sfun, ...)
  chain <- matrix(NA_real_, n_pilot, d)
  s_d <- (2.38^2) / d
  scale_factor <- 1
  shape_sd <- init_sd            # relative per-parameter shape
  n_accept_block <- 0L
  block_accepts <- numeric(0)

  for (it in 1:n_pilot) {
    step_sd <- sqrt(scale_factor * s_d) * shape_sd
    prop <- theta + rnorm(d, 0, step_sd)
    prp <- logpost(prop, dat, Sfun, ...)
    if (is.finite(prp$lp) && log(runif(1)) < (prp$lp - cur$lp)) {
      theta <- prop; cur <- prp; n_accept_block <- n_accept_block + 1L
    }
    chain[it, ] <- theta

    if (it %% adapt_every == 0) {
      ## WINDOWED, not cumulative -- only the most recent `window`
      ## draws inform the shape, so the initial transient (theta0 -> the
      ## posterior mode) can't contaminate the final per-parameter scale.
      lo <- max(1, it - window + 1)
      recent_sd <- apply(chain[lo:it, , drop = FALSE], 2, sd)
      recent_sd[!is.finite(recent_sd) | recent_sd < 1e-6] <- shape_sd[!is.finite(recent_sd) | recent_sd < 1e-6]
      shape_sd <- recent_sd

      block_accept <- n_accept_block / adapt_every
      block_accepts <- c(block_accepts, block_accept)
      scale_factor <- scale_factor * exp((block_accept - target_accept) * 0.5)
      scale_factor <- min(max(scale_factor, 1e-4), 1e4)
      n_accept_block <- 0L
    }
  }
  cat("Pilot block acceptance rates (last 5 blocks): ",
      paste(round(tail(block_accepts, 5), 3), collapse = ", "), "\n")
  final_step_sd <- sqrt(scale_factor * s_d) * shape_sd
  list(step_sd = final_step_sd, last_theta = theta, block_accepts = block_accepts)
}

## Short validation probe with the FROZEN step sizes -- if realized
## acceptance is far from target, rescale (uniformly across all
## coordinates, since we're not re-estimating shape here) and recheck.
calibrate_diagonal <- function(dat, Sfun, theta_start, step_sd,
                                n_probe = 1000, target_accept = 0.234,
                                tol = 0.07, max_rounds = 4, ...) {
  sd_vec <- step_sd
  theta <- theta_start
  for (round in seq_len(max_rounds)) {
    probe <- run_mh(dat, Sfun, n_iter = n_probe, theta0 = theta,
                     step_sd = sd_vec, verbose_every = 0, ...)
    cat(sprintf("Calibration round %d: accept_rate=%.3f\n", round, probe$accept_rate))
    theta <- probe$chain[n_probe, ]
    if (abs(probe$accept_rate - target_accept) < tol) break
    ratio <- max(probe$accept_rate, 1e-3) / target_accept
    sd_vec <- sd_vec * ratio^(2 / length(step_sd))
  }
  list(step_sd = sd_vec, theta_start = theta)
}

run_mh <- function(dat, Sfun, n_iter = 6000, theta0 = rep(0, 9),
                    step_sd, verbose_every = 0, ...) {
  d <- length(theta0)
  theta <- theta0
  cur <- logpost(theta, dat, Sfun, ...)
  chain <- matrix(NA_real_, n_iter, d)
  n_accept <- 0L
  total_terms <- 0
  for (it in 1:n_iter) {
    prop <- theta + rnorm(d, 0, step_sd)
    prp <- logpost(prop, dat, Sfun, ...)
    total_terms <- total_terms + prp$n_terms
    if (is.finite(prp$lp) && log(runif(1)) < (prp$lp - cur$lp)) {
      theta <- prop; cur <- prp; n_accept <- n_accept + 1L
    }
    chain[it, ] <- theta
    if (verbose_every > 0 && it %% verbose_every == 0) {
      cat(sprintf("  iter %d / %d, accept rate so far: %.2f\n",
                   it, n_iter, n_accept / it))
    }
  }
  list(chain = chain, accept_rate = n_accept / n_iter,
       mean_terms_per_eval = total_terms / (n_iter * nrow(dat)))
}

par_names <- c("alpha0","alpha1","alpha2","beta0","beta1","beta2",
               "gamma0","gamma1","gamma2")
theta0 <- c(0, 0, 0, log(mean(dat$roots[dat$roots > 0])), 0, 0, 0, 0, 0)

## ---------------------------------------------------------------------
## 4. Tune + calibrate the diagonal step sizes, then run all three
##    truncation strategies with the SAME step sizes.
## ---------------------------------------------------------------------
burn    <- 2000          # paper doesn't state burn-in length explicitly;
                          # kept generous given the ridge in alpha0/alpha2
n_post  <- 40000         # <- matches "40,000 MCMC posterior samples after burn-in"
thin    <- 20            # <- matches "every twentieth sample"
n_iter  <- burn + n_post
K_fixed <- 100
verbose_every <- 200
eps <- 2^-52 * 10^6

cat("\nTuning diagonal step sizes (windowed, no correlation estimation)...\n")
set.seed(0)
pilot <- tune_diagonal_scale(dat, S_fixed, theta0 = theta0, n_pilot = 3000, K = 100)

cat("\nCalibrating frozen step sizes against the target acceptance rate...\n")
set.seed(0)
calib <- calibrate_diagonal(dat, S_fixed, theta_start = pilot$last_theta,
                             step_sd = pilot$step_sd, K = 100)
step_sd <- calib$step_sd
theta_start <- calib$theta_start
cat("\nstep_sd = ", round(step_sd, 4), "\n")
cat("theta_start = ", round(theta_start, 4), "\n")

cat("\nRunning MH with fixed truncation (K =", K_fixed, ", matches the paper)...\n")
set.seed(1)
t_fixed <- system.time(
  fit_fixed <- run_mh(dat, S_fixed, n_iter = n_iter, theta0 = theta_start,
                       step_sd = step_sd, verbose_every = verbose_every, K = K_fixed)
)

cat("Running MH with Sum-to-threshold truncation (eps = 2.2e-10)...\n")
set.seed(1)
t_thresh <- system.time(
  fit_thresh <- run_mh(dat, S_threshold, n_iter = n_iter, theta0 = theta_start,
                        step_sd = step_sd, verbose_every = verbose_every, eps = eps)
)

cat("Running MH with Error-bounding pairs truncation (eps = 2.2e-10)...\n")
set.seed(1)
t_bound <- system.time(
  fit_bound <- run_mh(dat, S_bounding, n_iter = n_iter, theta0 = theta_start,
                       step_sd = step_sd, verbose_every = verbose_every, eps = eps)
)

## FIX: a second, independently-started/seeded chain (Fixed K=100 only,
## for speed) so we get Gelman-Rubin R-hat -- a much stronger convergence
## check than Geweke on a single chain.
cat("\nRunning a second, independently-started chain for R-hat...\n")
set.seed(99)
fit_fixed_chain2 <- run_mh(dat, S_fixed, n_iter = n_iter, theta0 = rep(0, 9),
                            step_sd = step_sd, verbose_every = 0, K = 100)

## ---------------------------------------------------------------------
## 5. Posterior summaries (thinned, matching the paper's approach)
## ---------------------------------------------------------------------
keep <- seq(burn + 1, n_iter, by = thin)

summarise_chain <- function(fit, name) {
  post <- fit$chain[keep, , drop = FALSE]
  data.frame(method = name, parameter = par_names,
             mean = colMeans(post), sd = apply(post, 2, sd))
}

results <- rbind(
  summarise_chain(fit_fixed,  "Fixed K=100"),
  summarise_chain(fit_thresh, "Sum-to-threshold"),
  summarise_chain(fit_bound,  "Error-bounding pairs")
)

cat("\n==================== Posterior summaries ====================\n")
print(reshape(results, idvar = "parameter", timevar = "method",
              direction = "wide"), row.names = FALSE)

cat("\n==================== Timing & acceptance ====================\n")
timing <- data.frame(
  method            = c("Fixed K=100", "Sum-to-threshold", "Error-bounding pairs"),
  time_sec          = c(t_fixed[3], t_thresh[3], t_bound[3]),
  accept_rate       = c(fit_fixed$accept_rate, fit_thresh$accept_rate, fit_bound$accept_rate),
  mean_terms_per_ev = c(fit_fixed$mean_terms_per_eval, fit_thresh$mean_terms_per_eval,
                         fit_bound$mean_terms_per_eval)
)
print(timing, row.names = FALSE)

## ---------------------------------------------------------------------
## 6a. Convergence diagnostics: Geweke (single chain) + Gelman-Rubin (2 chains)
## ---------------------------------------------------------------------
if (requireNamespace("coda", quietly = TRUE)) {
  mc1 <- coda::mcmc(fit_fixed$chain[keep, ])
  colnames(mc1) <- par_names
  cat("\n==================== Geweke diagnostic (Fixed K=100) ====================\n")
  print(coda::geweke.diag(mc1))

  mc2 <- coda::mcmc(fit_fixed_chain2$chain[keep, ])
  colnames(mc2) <- par_names
  cat("\n==================== Gelman-Rubin R-hat (2 independent chains) ====================\n")
  print(coda::gelman.diag(coda::mcmc.list(mc1, mc2), autoburnin = FALSE))
  cat("(R-hat should be close to 1.00, conventionally < 1.05, for every parameter.)\n")
} else {
  cat("\n(install.packages('coda') to get Geweke / R-hat diagnostics here)\n")
}

## ---------------------------------------------------------------------
## 6b. Effective Sample Size (ESS) and ESS/sec per truncation method
## ---------------------------------------------------------------------
if (requireNamespace("coda", quietly = TRUE)) {

  ess_table <- function(fit, time_sec, name) {
    post <- coda::mcmc(fit$chain[keep, , drop = FALSE])
    colnames(post) <- par_names
    ess <- coda::effectiveSize(post)
    data.frame(method      = name,
               parameter   = par_names,
               ESS         = round(ess, 1),
               ESS_per_sec = round(ess / time_sec, 2))
  }

  ess_all <- rbind(
    ess_table(fit_fixed,  t_fixed[3],  "Fixed K=100"),
    ess_table(fit_thresh, t_thresh[3], "Sum-to-threshold"),
    ess_table(fit_bound,  t_bound[3],  "Error-bounding pairs")
  )

  cat("\n==================== ESS and ESS/sec by parameter ====================\n")
  print(ess_all, row.names = FALSE)

  ## Worst-case (min) and typical (median) ESS/sec per method -- the min
  ## across parameters is what actually limits how long you need to run.
  ess_summary <- do.call(rbind, lapply(split(ess_all, ess_all$method), function(d) {
    data.frame(method          = d$method[1],
               min_ESS         = round(min(d$ESS), 1),
               median_ESS      = round(median(d$ESS), 1),
               min_ESS_per_sec = round(min(d$ESS_per_sec), 3),
               med_ESS_per_sec = round(median(d$ESS_per_sec), 3))
  }))

  cat("\n==================== ESS/sec summary (worst / typical parameter) ====================\n")
  print(ess_summary, row.names = FALSE)

} else {
  cat("\n(install.packages('coda') to compute ESS / ESS-per-second here)\n")
}

## ---------------------------------------------------------------------
## 7. Independent sanity check (NOT derived from the paper's Table 4):
##    does the fitted zero-inflation probability match the raw observed
##    zero rates by photoperiod?
## ---------------------------------------------------------------------
theta_hat <- colMeans(fit_fixed$chain[keep, , drop = FALSE])
p_8h  <- plogis(theta_hat[1])
p_16h <- plogis(theta_hat[1] + theta_hat[3])
cat("\n==================== Independent sanity check ====================\n")
cat(sprintf("Fitted structural-zero prob: 8h = %.4f (raw zero rate %.4f)\n", p_8h, raw_zero_8h))
cat(sprintf("Fitted structural-zero prob: 16h = %.4f (raw zero rate %.4f)\n", p_16h, raw_zero_16h))
cat("(These won't match exactly -- some zeros come from the COM-Poisson\n",
    " component itself, not the structural-zero mechanism -- but fitted\n",
    " p should be well below the raw rate, not above it, and in the same\n",
    " ballpark, especially at 8h where mu is large and the COM-Poisson's\n",
    " own P(Y=0) is negligible.)\n")

## ---------------------------------------------------------------------
## 8. Was K = 100 enough? (unchanged; already established this isn't the bug)
## ---------------------------------------------------------------------
post_bound <- fit_bound$chain[keep, , drop = FALSE]
theta_hat_b <- colMeans(post_bound)

cells <- expand.grid(bap = c(2.2, 4.4, 8.8, 17.6), photo01 = c(0, 1))
cells$mu <- with(cells, exp(theta_hat_b[4] + theta_hat_b[5]*bap + theta_hat_b[6]*photo01))
cells$nu <- with(cells, exp(theta_hat_b[7] + theta_hat_b[8]*bap + theta_hat_b[9]*photo01))
cells$S_fixed100 <- mapply(function(mu, nu) S_fixed(mu, nu, K = 100)$S, cells$mu, cells$nu)
cells$S_exact    <- mapply(function(mu, nu) S_bounding(mu, nu, eps = 2.2e-16)$S, cells$mu, cells$nu)
cells$n_needed   <- mapply(function(mu, nu) S_bounding(mu, nu, eps = 2.2e-16)$n_terms, cells$mu, cells$nu)
cells$rel_error  <- abs(cells$S_fixed100 - cells$S_exact) / cells$S_exact

cat("\n============ Was K = 100 enough at the fitted (mu, nu)? ============\n")
print(cells[, c("bap","photo01","mu","nu","n_needed","rel_error")], row.names = FALSE, digits = 4)

## ---------------------------------------------------------------------
## 9. Fitted frequency distribution (cf. Table 6)
## ---------------------------------------------------------------------
compute_fitted_freq <- function(theta, dat, Sfun, kmax = 14, ...) {
  a0<-theta[1]; a1<-theta[2]; a2<-theta[3]
  b0<-theta[4]; b1<-theta[5]; b2<-theta[6]
  g0<-theta[7]; g1<-theta[8]; g2<-theta[9]
  p   <- plogis(a0 + a1*dat$bap + a2*dat$photo01)
  mu  <- exp(b0 + b1*dat$bap + b2*dat$photo01)
  nu  <- exp(g0 + g1*dat$bap + g2*dat$photo01)
  n <- nrow(dat)
  freq <- numeric(kmax + 2)
  for (i in seq_len(n)) {
    S <- Sfun(mu[i], nu[i], ...)$S
    p0 <- p[i] + (1 - p[i]) / S
    freq[1] <- freq[1] + p0
    cum <- p0
    for (k in 1:kmax) {
      pk <- (1 - p[i]) * exp(nu[i] * (k*log(mu[i]) - lgamma(k + 1))) / S
      freq[k + 1] <- freq[k + 1] + pk
      cum <- cum + pk
    }
    freq[kmax + 2] <- freq[kmax + 2] + max(0, 1 - cum)
  }
  freq
}

theta_fixed  <- theta_hat
theta_thresh <- colMeans(fit_thresh$chain[keep, , drop = FALSE])
theta_bound  <- theta_hat_b

freq_fixed  <- compute_fitted_freq(theta_fixed,  dat, S_fixed,     K = 100)
freq_thresh <- compute_fitted_freq(theta_thresh, dat, S_threshold, eps = 2.2e-10)
freq_bound  <- compute_fitted_freq(theta_bound,  dat, S_bounding,  eps = 2.2e-10)

observed <- c(64,10,13,15,21,18,24,21,23,21,17,12,5,2,3,1)

tab6 <- data.frame(
  Count = c(0:14, ">=15"), Observed = observed,
  `Fixed K=100` = round(freq_fixed, 1),
  `Sum-to-threshold` = round(freq_thresh, 1),
  `Error-bounding pairs` = round(freq_bound, 1),
  check.names = FALSE
)
cat("\n==================== Fitted frequency distribution (cf. Table 6) ====================\n")
print(tab6, row.names = FALSE)

chisq_stat <- function(obs, fit) sum((obs - fit)^2 / fit)
df_chisq <- length(observed) - 1 - length(theta_fixed)
chisq_tab <- data.frame(
  method = c("Fixed K=100", "Sum-to-threshold", "Error-bounding pairs"),
  chi2 = c(chisq_stat(observed, freq_fixed), chisq_stat(observed, freq_thresh),
           chisq_stat(observed, freq_bound)),
  df = df_chisq
)
chisq_tab$p_value <- pchisq(chisq_tab$chi2, chisq_tab$df, lower.tail = FALSE)
cat("\n")
print(chisq_tab, row.names = FALSE, digits = 4)