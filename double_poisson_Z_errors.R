  #!/usr/bin/env Rscript
  library(cmdstanr)
  library(future.apply)
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
    #if (k %% 1000 == 0) cat("K=", k, " ")
    
    if (k == 0) {
        1L
    } else {
        k <- as.numeric(k)
        
        -k +
        k * log(k) -
        lgamma(k + 1L) +
        theta[[2L]] * k +
        theta[[2L]] * k * log(theta[[1L]]) -
        theta[[2L]] * k * log(k)
    }
}

  format_result <- function(x) {
    if (inherits(x, "mpfr")) {
      x <- Rmpfr::asNumeric(x)
    }
    if (is.infinite(x) && x < 0) return("-Inf")
    formatC(x, format = "g", digits = 8)
  }

  # Use multicore (fork) on Unix/Mac if available, otherwise multisession
  if (.Platform$OS.type == "unix") {
    plan(multicore, workers = min(4, availableCores()))
  } else {
    plan(multisession, workers = min(4, availableCores()))
  }

  run_comparison <- function() {
    mu <- c(10, 100, 1000, 10000)
    nu <- c(0.1, 0.01, 0.001, 0.0001)
    lambda <- mu^nu
    log_lambda <- log(lambda)
    M <- c(10^5L, 10^5L, 10^5L, 10^5L)
    initial_k <- 0L
    prec <- 128L
    eps <- 2^-52

    # Each worker computes fixed_value ONCE and reuses it for both error
    # thresholds (2.2x10^-10 and 2.2x10^-16), mirroring the Python but
    # avoiding the redundant recomputation the Python does per block.
    errors <- future_lapply(seq_along(mu), function(i) {
      fixed_value <- fixed(f, c(log_lambda[i], nu[i]), M[i], initial_k,
                            prec = prec)$value

      seq_value_10 <- sequential(f, c(log_lambda[i], nu[i]), M[i],
                                  eps * 10^6, initial_k, prec = prec)$value
      bp_value_10 <- bounding_pairs(f, c(log_lambda[i], nu[i]), M[i], 0,
                                    eps * 10^6, initial_k, bucket_size = 1L,
                                    prec = prec)$value

      seq_value_16 <- sequential(f, c(log_lambda[i], nu[i]), M[i],
                                  eps, initial_k, prec = prec)$value
      bp_value_16 <- bounding_pairs(f, c(log_lambda[i], nu[i]), M[i], 0,
                                    eps, initial_k, bucket_size = 1L,
                                    prec = prec)$value

      list(
        fixed_value = fixed_value,
        m10 = list(
          sequential = exp(logdiffexp(seq_value_10, fixed_value, prec=prec)),
          BP         = exp(logdiffexp(bp_value_10, fixed_value, prec=prec))
        ),
        m16 = list(
          sequential = exp(logdiffexp(seq_value_16, fixed_value, prec=prec)),
          BP         = exp(logdiffexp(bp_value_16, fixed_value, prec=prec))
        )
      )
    }, future.seed = TRUE, future.packages = c("Rmpfr"))

    fixed_values   <- lapply(errors, `[[`, "fixed_value")
    error_minus_10 <- lapply(errors, `[[`, "m10")
    error_minus_16 <- lapply(errors, `[[`, "m16")

    rows <- lapply(seq_along(mu), function(i) {
      c(sprintf("mu=%g | phi=%g", mu[i], nu[i]),
        format_result(error_minus_10[[i]][["sequential"]]),
        format_result(error_minus_10[[i]][["BP"]]),
        format_result(error_minus_16[[i]][["sequential"]]),
        format_result(error_minus_16[[i]][["BP"]]))
    })

    headers <- c("", "2.2x10^-10|Sequential", "2.2x10^-10|BP",
                "2.2x10^-16|Sequential", "2.2x10^-16|BP")
    table <- do.call(rbind, rows)
    colnames(table) <- headers
    print(noquote(table), quote = FALSE, right = TRUE)
    invisible(table)
  }

  if (sys.nframe() == 0L) {
    run_comparison()
  }