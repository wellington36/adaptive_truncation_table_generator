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
  stop("Please run install.packages('agridat') first (CRAN), then ",
       "re-run this script, or set data_source <- 'table2'.")
}
data(ridout.appleshoots, package = "agridat")
dat <- ridout.appleshoots
dat$photo01 <- as.numeric(dat$photo == 16)  # 0 = 8h, 1 = 16h (matches x2 in the paper)


cat(sprintf("n = %d, %% zero = %.1f%%\n", nrow(dat), 100 * mean(dat$roots == 0)))

## ---------------------------------------------------------------------
## 1. Truncation strategies for S(mu, nu) = sum_s (mu^s/s!)^nu
##
## All three use the numerically stable ratio
##   r(n) = a_n / a_{n-1} = mu / n^nu
## instead of forming a_n directly, since this avoids overflow for
## large mu/small nu and never needs to divide two possibly-huge numbers.
## ---------------------------------------------------------------------

## (a) Fixed truncation -- what the original paper used (K = 100, i.e.
##     summing the first 101 terms, "truncated at the 101st term").
S_fixed <- function(mu, nu, K = 100) {
  ## a_n / a_{n-1} = (mu/n)^nu  -- NOTE: the whole ratio (mu/n) is raised
  ## to the nu power (this is Barriga & Louzada's mu-reparametrised
  ## COM-Poisson, Eq. (3). Worked in log space to stay stable across
  ## the huge nu/mu ranges MCMC proposals can transiently visit.
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
  nu <- exp(g0 + g1 * dat$bap + g2 * dat$photo01)

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
## 3. Random-walk Metropolis-Hastings sampler
##    (mirrors the paper's MH design: fixed proposal covariance, run
##    until convergence; here a diagonal random walk for simplicity)
## ---------------------------------------------------------------------
run_mh <- function(dat, Sfun, n_iter = 6000, prop_sd = default_prop_sd,
                    theta0 = rep(0, 9), verbose_every = 0, ...) {
  theta <- theta0
  cur <- logpost(theta, dat, Sfun, ...)
  chain <- matrix(NA_real_, n_iter, 9)
  n_accept <- 0L
  total_terms <- 0
  for (it in 1:n_iter) {
    prop <- theta + rnorm(9, 0, prop_sd)
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

## Per-parameter random-walk step sizes. A single shared step size mixes
## poorly here because the parameters live on very different scales (e.g.
## beta1's posterior SD is ~0.006 in the paper's Table 4, vs ~0.5 for
## alpha0/alpha2) -- with one shared prop_sd you either barely move the
## small-scale parameters or blow the acceptance rate on the large-scale
## ones. These are rough guesses based on Table 4 of the paper, scaled
## down a bit to target the ~15-40% acceptance range typical for
## random-walk MH
default_prop_sd <- c(alpha0=0.25, alpha1=0.02, alpha2=0.25,
                      beta0=0.04, beta1=0.004, beta2=0.12,
                      gamma0=0.12, gamma1=0.012, gamma2=0.18)

## Starting values: rough, data-driven guess (feel free to replace with
## Table 4 of the paper if you want faster burn-in on the real data:
## alpha=(3.57,-0.002,3.30) [sign convention differs -- see note below],
## beta=(1.85,0.011,-0.65), gamma=(-0.80,0.059,-0.90)).
theta0 <- c(0, 0, 0, log(mean(dat$roots[dat$roots > 0])), 0, 0, 0, 0, 0)

## ---------------------------------------------------------------------
## 4. Run the same sampler under the three truncation strategies
## ---------------------------------------------------------------------
n_iter <- 6000   # paper used 40,000 (thinned); raise this once you've
                 # confirmed the pipeline runs cleanly on your machine
burn   <- 1500

cat("\nRunning MH with fixed truncation (K = 100, matches the paper)...\n")
set.seed(1)
t_fixed <- system.time(
  fit_fixed <- run_mh(dat, S_fixed, n_iter = n_iter, theta0 = theta0, verbose_every = 500, K = 100)
)

cat("Running MH with Sum-to-threshold truncation (eps = 2.2e-10)...\n")
set.seed(1)
t_thresh <- system.time(
  fit_thresh <- run_mh(dat, S_threshold, n_iter = n_iter, theta0 = theta0, verbose_every = 500, eps = 2.2e-10)
)

cat("Running MH with Error-bounding pairs truncation (eps = 2.2e-10)...\n")
set.seed(1)
t_bound <- system.time(
  fit_bound <- run_mh(dat, S_bounding, n_iter = n_iter, theta0 = theta0, verbose_every = 500, eps = 2.2e-10)
)

## ---------------------------------------------------------------------
## 5. Compare posterior summaries, timing and truncation cost
## ---------------------------------------------------------------------
summarise_chain <- function(fit, name) {
  post <- fit$chain[-(1:burn), , drop = FALSE]
  data.frame(
    method    = name,
    parameter = par_names,
    mean      = colMeans(post),
    sd        = apply(post, 2, sd)
  )
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
## 6. Direct diagnostic: was K = 100 actually enough?
## Evaluate S(mu, nu) at the posterior-mean (mu_i, nu_i) for each of
## the 8 BAP x photoperiod treatment cells, comparing the paper's fixed
## K = 100 truncation against the guaranteed-accurate Error-bounding
## pairs result at eps = 2.2e-16.
## ---------------------------------------------------------------------
post_bound <- fit_bound$chain[-(1:burn), , drop = FALSE]
theta_hat <- colMeans(post_bound)

cells <- expand.grid(bap = c(2.2, 4.4, 8.8, 17.6), photo01 = c(0, 1))
cells$mu  <- with(cells, exp(theta_hat[4] + theta_hat[5]*bap + theta_hat[6]*photo01))
cells$nu <- with(cells, exp(theta_hat[7] + theta_hat[8]*bap + theta_hat[9]*photo01))

cells$S_fixed100  <- mapply(function(mu, nu) S_fixed(mu, nu, K = 100)$S, cells$mu, cells$nu)
cells$S_exact     <- mapply(function(mu, nu) S_bounding(mu, nu, eps = 2.2e-16)$S, cells$mu, cells$nu)
cells$n_needed    <- mapply(function(mu, nu) S_bounding(mu, nu, eps = 2.2e-16)$n_terms, cells$mu, cells$nu)
cells$rel_error   <- abs(cells$S_fixed100 - cells$S_exact) / cells$S_exact

cat("\n============ Was K = 100 enough at the fitted (mu, nu)? ============\n")
print(cells[, c("bap","photo01","mu","nu","n_needed","rel_error")], row.names = FALSE,
      digits = 4)
cat("\n(n_needed = terms required by Error-bounding pairs for eps = 2.2e-16;\n",
    " rel_error = relative error of the paper's K = 100 truncation at that cell.)\n")

## ---------------------------------------------------------------------
## 7. Fitted frequency distribution, analogous to Table 6 of Barriga &
##    Louzada (2014), but comparing the three truncation strategies
##    against each other (rather than ZICOM-Poisson vs. ZIP vs. ZINB as
##    the original table does).
##
## For each individual i, the fitted probability of Y_i = k is computed
## under the ZICOM-Poisson pmf (Eq. zicmp_pmf) at that method's posterior
## mean, using that method's own truncation of S(mu_i, nu_i); expected
## frequency for count k is the sum of these probabilities over all n
## individuals, following the population-averaged approach used to build
## Table 6. Counts 15+ are pooled into a single ">=15" bin, matching the
## paper's table.
## ---------------------------------------------------------------------
compute_fitted_freq <- function(theta, dat, Sfun, kmax = 14, ...) {
  a0<-theta[1]; a1<-theta[2]; a2<-theta[3]
  b0<-theta[4]; b1<-theta[5]; b2<-theta[6]
  g0<-theta[7]; g1<-theta[8]; g2<-theta[9]
  p   <- plogis(a0 + a1*dat$bap + a2*dat$photo01)
  mu  <- exp(b0 + b1*dat$bap + b2*dat$photo01)
  nu <- exp(g0 + g1*dat$bap + g2*dat$photo01)
  n <- nrow(dat)
  freq <- numeric(kmax + 2)  # bins: 0, 1, ..., kmax, >=kmax+1
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

theta_fixed  <- colMeans(fit_fixed$chain[-(1:burn), , drop = FALSE])
theta_thresh <- colMeans(fit_thresh$chain[-(1:burn), , drop = FALSE])
theta_bound  <- colMeans(fit_bound$chain[-(1:burn), , drop = FALSE])

freq_fixed  <- compute_fitted_freq(theta_fixed,  dat, S_fixed,     K = 100)
freq_thresh <- compute_fitted_freq(theta_thresh, dat, S_threshold, eps = 2.2e-10)
freq_bound  <- compute_fitted_freq(theta_bound,  dat, S_bounding,  eps = 2.2e-10)

## Observed frequency distribution, taken directly from Table 6 of
## Barriga & Louzada (2014) (independent of our own Table 2
## reconstruction, so unaffected by the small n=269-vs-270 discrepancy
## noted in Section 0) -- REPLACE with your own tabulated `dat$roots`
## frequencies if you are working from verified raw data instead.
observed <- c(64,10,13,15,21,18,24,21,23,21,17,12,5,2,3,1)

tab6 <- data.frame(
  Count               = c(0:14, ">=15"),
  Observed            = observed,
  `Fixed K=100`        = round(freq_fixed, 1),
  `Sum-to-threshold`   = round(freq_thresh, 1),
  `Error-bounding pairs` = round(freq_bound, 1),
  check.names = FALSE
)

cat("\n==================== Fitted frequency distribution (cf. Table 6) ====================\n")
print(tab6, row.names = FALSE)

chisq_stat <- function(obs, fit) sum((obs - fit)^2 / fit)
df_chisq <- length(observed) - 1 - length(theta_fixed)  # 16 categories - 1 - 9 regr. params

chisq_tab <- data.frame(
  method  = c("Fixed K=100", "Sum-to-threshold", "Error-bounding pairs"),
  chi2    = c(chisq_stat(observed, freq_fixed),
              chisq_stat(observed, freq_thresh),
              chisq_stat(observed, freq_bound)),
  df      = df_chisq
)
chisq_tab$p_value <- pchisq(chisq_tab$chi2, chisq_tab$df, lower.tail = FALSE)

cat("\n")
print(chisq_tab, row.names = FALSE, digits = 4)
cat("\n(df = 16 categories - 1 - 9 regression parameters; the paper's own df\n",
    " convention for Table 6 is not stated, so treat this as a self-consistent\n",
    " comparison ACROSS our three truncation methods rather than a like-for-like\n",
    " replication of the paper's own chi2/p-value column.)\n")