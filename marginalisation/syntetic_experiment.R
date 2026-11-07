suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(progress)
})

# auxiliary functions
# log(exp(x) + exp(y)), -Inf safe
logspace_add <- function(x, y) {
  m <- pmax(x, y)
  out <- m + log1p(exp(-abs(x - y)))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}

log_sum_exp <- function(x) {
  m <- max(x); if (is.infinite(m) && m < 0) return(-Inf)
  m + log(sum(exp(x - m)))
}

# log(1 - exp(x)), x <= 0
log1mexp <- function(x) {
  if (any(x > 0, na.rm = TRUE)) stop("log1mexp requires x <= 0")
  ifelse(x > -log(2), log(-expm1(x)), log1p(-exp(x)))
}

log1pexp <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

# log(log(1 + exp(x))), no underflow
log_log1pexp <- function(x) {
  out <- x
  big <- x >= -30
  out[big] <- log(log1pexp(x[big]))
  out
}

# log PMF of total progeny, NB offspring
log_rj <- function(j, R0, k) {
  lgamma((k + 1) * j - 1) - lgamma(k * j) - lgamma(j + 1) +
    (j - 1) * (log(R0) - log(k)) - ((k + 1) * j - 1) * log1p(R0 / k)
}

# log of the limiting ratio L = q * R0 * ((k+1)/(k+R0))^(k+1); (k+1) = (k+R0) + (1-R0)
log_L_fun <- function(R0, k, log_q) log_q + log(R0) + (k + 1) * log1p((1 - R0) / (k + R0))

# Adaptive enclosure of rho_0 (all arithmetic in log scale)
# Uses a_n = t_{n+1}/t_n (Algorithm 1): lower tail t_n a_n/(1-a_n), upper t_n L/(1-L).
# Returns log(1 - rho) at both endpoints, never rho itself.
adaptive_rho <- function(R0, k, p, eps, m, max_n = 2e5L, b0 = 64L, b_max = 4096L) {
  log_q  <- log1p(-p)
  log_L  <- log_L_fun(R0, k, log_q)
  log1mL <- log1mexp(log_L)
  log_eps <- log(eps)

  n <- 0L; log_S <- -Inf; b <- b0
  repeat {
    js  <- (n + 1L):(n + b)
    lr  <- log_rj(n + 1L + 0:b, R0, k)             # one extra term for the ratio
    log_t <- lr[1:b] + js * log_q
    log_a <- log_q + lr[2:(b + 1L)] - lr[1:b]      # log a_j

    # partial sums S_j for the whole block (cumsum uses long double in R)
    mx <- max(log_t)
    log_S_c <- logspace_add(log_S, mx + log(cumsum(exp(log_t - mx))))

    ok <- is.finite(log_a) & (log_a < log_L)       # ratio resolved below its limit
    log_w <- rep(Inf, b)
    if (any(ok)) {
      la <- log_a[ok]; lt <- log_t[ok]; lS <- log_S_c[ok]
      log1ma   <- log1mexp(la)
      log_gap  <- log_L + log1mexp(la - log_L)                    # log(L - a)
      log_lo   <- logspace_add(lS, lt + la - log1ma)
      log_hi   <- logspace_add(lS, lt + log_L - log1mL)
      # hi - lo = t (L - a) / ((1-L)(1-a)), computed without subtracting endpoints
      log_D    <- lt + log_gap - log1mL - log1ma
      good <- log_hi < 0
      lw <- rep(Inf, length(la))
      if (any(good)) {
        l1m_hi <- log1mexp(log_hi[good])
        lw[good] <- log(m) + log_log1pexp(log_D[good] - l1m_hi)   # log of m*log1p(D/(1-hi))
      }
      log_w[ok] <- lw
      hit <- which(log_w <= log_eps)
      if (length(hit)) {
        i  <- hit[1]; idx <- which(ok)[ ok_index <- which(lw <= log_eps)[1] ]
        l1m_lo <- log1mexp(log_lo[ok_index]); l1m_hi <- log1mexp(log_hi[ok_index])
        return(list(n = js[idx], log1m_lo = l1m_lo, log1m_hi = l1m_hi,
                    log1m_mid = 0.5 * (l1m_lo + l1m_hi), width = exp(lw[ok_index])))
      }
    }
    n <- n + b; log_S <- log_S_c[b]
    if (n >= max_n) stop("max_n reached before tolerance")
    b <- min(2L * b, b_max)
  }
}

# Brute-force partial sum, J terms, log scale. Returns log(1 - rho_J).
log1m_rho_fixed <- function(R0, k, p, J) {
  j <- seq_len(J)
  lrho <- log_sum_exp(log_rj(j, R0, k) + j * log1p(-p))
  if (lrho >= 0) -Inf else log1mexp(lrho)
}

# Data + likelihood (sufficient statistics: unique sizes and counts)
make_data <- function(y) {
  tb <- table(y); list(u = as.numeric(names(tb)), cnt = as.numeric(tb), m = length(y))
}
ll_const <- function(d, R0, k, p)
  sum(d$cnt * (log_rj(d$u, R0, k) + log1mexp(d$u * log1p(-p))))

loglik_at <- function(d, R0, k, p, type, eps = NA, J = NA, max_n = 2e5L) {
  cst <- ll_const(d, R0, k, p)
  if (type == "adaptive") {
    z <- adaptive_rho(R0, k, p, eps, d$m, max_n)
    list(ll = cst - d$m * z$log1m_mid,
         ll_lo = cst - d$m * z$log1m_lo,            # smaller rho -> smaller ll
         ll_hi = cst - d$m * z$log1m_hi,
         terms = z$n, width = z$width)
  } else {
    l1m <- log1m_rho_fixed(R0, k, p, J)
    list(ll = cst - d$m * l1m, ll_lo = NA_real_, ll_hi = NA_real_,
         terms = J, width = NA_real_)
  }
}

# Detection-conditioned sampler; table length chosen with a rigorous tail bound
# (r_{j+1}/r_j increases to f <= 1, so tail_{>J} <= r_J f/(1-f))
make_sampler <- function(R0, k, p, tail_tol = 1e-14, J = 2000L, J_max = 1e7L) {
  log_f <- log(R0) + (k + 1) * log1p((1 - R0) / (k + R0))
  repeat {
    j  <- seq_len(J)
    lw <- log_rj(j, R0, k) + log1mexp(j * log1p(-p))
    lZ <- log_sum_exp(lw)
    log_tail <- log_rj(J, R0, k) + log_f - log1mexp(log_f)
    if (log_tail - lZ < log(tail_tol)) break
    J <- J * 2L; if (J > J_max) stop("sampler table too long")
  }
  list(prob = exp(lw - lZ), J = J)
}

# MLE with a given likelihood evaluator
BND <- list(R0_max = 0.99, lk = log(c(0.05, 50)), p_lo = 1e-3)
decode <- function(z) c(R0 = BND$R0_max * plogis(z[1]),
                        k  = exp(BND$lk[1] + diff(BND$lk) * plogis(z[2])),
                        p  = BND$p_lo + (1 - 2 * BND$p_lo) * plogis(z[3]))
encode <- function(th) c(qlogis(th[1] / BND$R0_max),
                         qlogis((log(th[2]) - BND$lk[1]) / diff(BND$lk)),
                         qlogis((th[3] - BND$p_lo) / (1 - 2 * BND$p_lo)))

fit_one <- function(d, spec, start = c(.5, 1, .2), max_n = 2e5L) {
  cnt <- new.env(); cnt$terms <- 0; cnt$evals <- 0
  obj <- function(z) {
    th <- decode(z)
    r <- tryCatch(loglik_at(d, th[1], th[2], th[3], spec$type, spec$eps, spec$J, max_n),
                  error = function(e) NULL)
    cnt$evals <- cnt$evals + 1
    if (is.null(r) || !is.finite(r$ll)) return(1e10)
    cnt$terms <- cnt$terms + r$terms
    -r$ll
  }
  t0 <- proc.time()[["elapsed"]]
  o <- optim(encode(start), obj, method = "Nelder-Mead",
             control = list(maxit = 3000, reltol = 1e-12))
  o <- optim(o$par, obj, method = "Nelder-Mead",
             control = list(maxit = 3000, reltol = 1e-12))
  el <- proc.time()[["elapsed"]] - t0
  th <- decode(o$par)
  fin <- tryCatch(loglik_at(d, th[1], th[2], th[3], spec$type, spec$eps, spec$J, max_n),
                  error = function(e) list(ll = NA_real_, ll_lo = NA_real_, ll_hi = NA_real_,
                                           terms = NA_real_, width = NA_real_))
  list(R0 = th[[1]], k = th[[2]], p = th[[3]], ll_own = fin$ll,
       ll_lo = fin$ll_lo, ll_hi = fin$ll_hi, terms_final = fin$terms, width = fin$width,
       terms_total = cnt$terms, evals = cnt$evals, conv = o$convergence, time = el)
}

# Syntetic experiment
run_experiment1 <- function(
    grid = expand_grid(R0 = c(.5, .7, .9), k = c(.1, 1, 10), p_sent = c(.01, .1, .5)),
    B = 50L, m = 500L, seed = 42L,
    adapt_eps = c(1e-2, 1e-4, 1e-8, 1e-10), fixed_J = c(50L, 200L, 1000L),
    eps_ref = 1e-12, slack = 1e-11, max_n = 2e5L, max_n_ref = 1e6L) {

  specs <- c(
    lapply(adapt_eps, function(e) list(name = sprintf("adaptive_%.0e", e),
                                       type = "adaptive", eps = e, J = NA)),
    lapply(fixed_J,   function(J) list(name = paste0("fixed_", J),
                                       type = "fixed", eps = NA, J = J)))
  ref_spec <- list(type = "adaptive", eps = eps_ref, J = NA)

  pb <- progress_bar$new(
    format = "  [:bar] :percent | :current/:total reps | elapsed :elapsed | eta :eta",
    total = nrow(grid) * B, clear = FALSE, width = 80)

  out <- list()
  for (g in seq_len(nrow(grid))) {
    P <- grid[g, ]
    samp <- make_sampler(P$R0, P$k, P$p_sent)

    for (b in seq_len(B)) {
      set.seed(seed + 1000L * g + b)
      d <- make_data(sample.int(samp$J, m, replace = TRUE, prob = samp$prob))

      ref <- fit_one(d, ref_spec, max_n = max_n_ref)

      # independent check of the reference: brute force with 3x the terms it used
      cst <- ll_const(d, ref$R0, ref$k, ref$p)
      ll_brute <- cst - m * log1m_rho_fixed(ref$R0, ref$k, ref$p, 3L * ref$terms_final)
      ref_brute_diff <- ref$ll_own - ll_brute
      if (!is.finite(ref_brute_diff) || abs(ref_brute_diff) > 10 * eps_ref + slack)
        pb$message(sprintf("REFERENCE MISMATCH cell=%d rep=%d R0=%.3f k=%.3f p=%.3f diff=%.3e",
                           g, b, P$R0, P$k, P$p_sent, ref_brute_diff))

      for (s in specs) {
        f <- fit_one(d, s, max_n = max_n)

        # reference likelihood evaluated at THIS method's estimate
        ll_ref_at <- tryCatch(
          loglik_at(d, f$R0, f$k, f$p, "adaptive", eps_ref, NA, max_n_ref)$ll,
          error = function(e) NA_real_)

        contains <- if (s$type == "adaptive")
          (ll_ref_at >= f$ll_lo - slack) && (ll_ref_at <= f$ll_hi + slack) else NA
        if (isFALSE(contains))
          pb$message(sprintf(
            "COUNTER-EXAMPLE %s cell=%d rep=%d R0=%.4f k=%.4f p=%.4f: [%.12f, %.12f] vs ref %.12f",
            s$name, g, b, f$R0, f$k, f$p, f$ll_lo, f$ll_hi, ll_ref_at))

        out[[length(out) + 1L]] <- tibble(
          R0 = P$R0, k = P$k, p_sent = P$p_sent, rep = b, method = s$name,
          est_R0 = f$R0, est_k = f$k, est_p = f$p,
          ref_R0 = ref$R0, ref_k = ref$k, ref_p = ref$p,
          ll_eval_err = f$ll_own - ll_ref_at,         # error of the likelihood at the estimate
          regret      = ref$ll_own - ll_ref_at,       # loss in true (reference) likelihood
          contains_ref = contains,
          width = f$width, terms_final = f$terms_final,
          terms_total = f$terms_total, evals = f$evals,
          time = f$time, conv = f$conv, ref_brute_diff = ref_brute_diff)
      }
      pb$tick()
    }
  }
  bind_rows(out)
}

summarise_experiment1 <- function(res) {
  res %>%
    group_by(R0, k, p_sent, method) %>%
    summarise(
      n = n(),
      mean_abs_ll_err = mean(abs(ll_eval_err), na.rm = TRUE),
      p95_abs_ll_err  = quantile(abs(ll_eval_err), .95, na.rm = TRUE),
      mean_regret     = mean(regret, na.rm = TRUE),
      p95_regret      = quantile(regret, .95, na.rm = TRUE),
      rmse_R0_vs_ref  = sqrt(mean((est_R0 - ref_R0)^2, na.rm = TRUE)),
      rmse_k_vs_ref   = sqrt(mean((est_k - ref_k)^2, na.rm = TRUE)),
      rmse_p_vs_ref   = sqrt(mean((est_p - ref_p)^2, na.rm = TRUE)),
      coverage        = mean(contains_ref, na.rm = TRUE),   # should be 1 for adaptive
      median_terms_total = median(terms_total, na.rm = TRUE),
      median_time     = median(time, na.rm = TRUE),
      .groups = "drop")
}

# ---------------- run ----------------
res <- run_experiment1(
  grid = expand_grid(R0 = c(.5, .7, .9), k = c(.1, 1, 10), p_sent = c(.01, .1, .5)),
  B = 100L, m = 500L,
  adapt_eps = c(1e-4, 1e-7, 1e-10),
  fixed_J   = c(50L, 200L, 1000L))

print(table(res$method, res$conv))   # 0 = optim converged, 1 = maxit reached
