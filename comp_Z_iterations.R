#!/usr/bin/env Rscript
library(cmdstanr)
library(COMPoissonReg)

if (interactive()) {
    # Running interactively in Positron/RStudio
    root <- getwd()
} else {
    # Running as Rscript from terminal
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)

    if (length(file_arg) > 0) {
        script_path <- sub("^--file=", "", file_arg[[1L]])
        root <- dirname(normalizePath(script_path, mustWork = FALSE))
    } else {
        root <- getwd()
    }
}

project_root <- root

source(file.path(project_root, "utils", "utils.R"))
source(file.path(project_root, "src", "fixed.R"))
source(file.path(project_root, "src", "sequential.R"))
source(file.path(project_root, "src", "bounding_pairs.R"))

f <- function(theta, k) {
  as.numeric(k) * theta[[1L]] - theta[[2L]] * lgamma(as.numeric(k) + 1)
}

format_result <- function(x) {
  if (is.infinite(x) && x < 0) return("-Inf")
  formatC(x, format = "g", digits = 8)
}

run_comparison <- function() {
  #mu <- c(10, 100, 1000, 10000)
  mu <- c(10, 100, 1000)
  #nu <- c(0.1, 0.01, 0.001, 0.0001)
  nu <- c(0.1, 0.01, 0.001)
  lambda <- mu^nu
  log_lambda <- log(lambda)
  #M <- c(10^4L, 10^5L, 10^5L, 3 * 10^5L)
  M <- c(10000L, 10000L, 10000L)
  initial_k <- 0L
  prec <- 128L
  eps <- 2^-10

  error_minus_10 <- lapply(seq_along(mu), function(i) {
    bp <- bounding_pairs(f, c(log_lambda[i], nu[i]), M[i], 0, eps * 10^6,
                         initial_k, bucket_size = 1L, prec=128L)$iterations
    sequential <- sequential(f, c(log_lambda[i], nu[i]), M[i],
                                  eps * 10^6, initial_k, prec=128L)$iterations
    c(sequential = sequential, BP = bp)
  })

  error_minus_16 <- lapply(seq_along(mu), function(i) {
    bp <- bounding_pairs(f, c(log_lambda[i], nu[i]), M[i], 0, eps,
                         initial_k, bucket_size = 1L, prec=128L)$iterations
    sequential <- sequential(f, c(log_lambda[i], nu[i]), M[i],
                                  eps, initial_k, prec=128L)$iterations
    c(sequential = sequential, BP = bp)
  })


  mod <- cmdstan_model(
    "stan/comp_Z_brms_fixed.stan",
    compile_standalone = TRUE,
    force_recompile = TRUE
  )
  mod$expose_functions(global = TRUE)

  libraries <- lapply(seq_along(mu), function(i) {
    fixed_value <- fixed(f, c(log_lambda[i], nu[i]), M[i], initial_k, prec = prec)$value

    brms_value <- log_z_com_poisson(log(mu[i]), nu[i])

    compoissonreg_log_pmf0 <- log(dcmp(0, lambda[i], nu[i]))
    
    c(brms = exp(logdiffexp(fixed_value, brms_value)),
      COMPoissonReg = exp(logdiffexp(fixed_value, -compoissonreg_log_pmf0)))
  })

  rows <- lapply(seq_along(mu), function(i) {
    c(sprintf("mu=%g | nu=%g", mu[i], nu[i]),
      error_minus_10[[i]][["sequential"]], error_minus_10[[i]][["BP"]],
      error_minus_16[[i]][["sequential"]], error_minus_16[[i]][["BP"]],
      print(libraries)
      format_result(libraries[[i]][["brms"]]),
      format_result(libraries[[i]][["COMPoissonReg"]]))
  })

  headers <- c("", "2.2x10^-10|Sequential", "2.2x10^-10|BP",
               "2.2x10^-16|Sequential", "2.2x10^-16|BP", "brms",
               "COMPoissonReg")
  table <- do.call(rbind, rows)
  colnames(table) <- headers
  print(noquote(table), quote = FALSE, right = TRUE)
  invisible(table)
}

if (sys.nframe() == 0L) {
  run_comparison()
}
