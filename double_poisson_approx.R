library(ggplot2)
library(gridExtra)

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
source(file.path(project_root, "src", "bounding_pairs.R"))
 
# ---- Double Poisson (Efron 1986) log unnormalised term ---------------------
# theta = c(mu, phi); matches gamlss.dist::DPO() parameterisation, i.e.
# p(y|mu,phi) = exp(-y) y^y / y! * (e*mu/y)^(phi*y) * sqrt(phi) * exp(-phi*mu)
log_term_dp <- function(theta, n) {
  mu <- theta[1]; phi <- theta[2]
  if (n == 0) return(0.5 * log(phi) - phi * mu)
  n <- as.numeric(n)
  n * log(n) - n - lgamma(n + 1) + phi * n * (1 + log(mu) - log(n)) + 0.5 * log(phi) - phi * mu
}
 
# L = 0 for the Double Poisson normalising constant (Appendix F derivation):
# a_{n+1}/a_n -> 0 as n -> infinity for all mu, phi > 0.
compute_K <- function(mu, phi, eps = 2.2e-16, M = 5000L, prec = 64L) {
  out <- bounding_pairs(f = log_term_dp, theta = c(mu, phi), M = M, L = 0,
                         eps = eps, initial_k = 0L, bucket_size = 1L, prec = prec)
  list(K = as.numeric(exp(out$value)), iterations = out$iterations)
}
 
# ---- grid -------------------------------------------------------------------
mus <- c(2, 5, 10)
phi_grid <- exp(seq(log(0.05), log(12), length.out = 80))
 
results <- do.call(rbind, lapply(mus, function(mu) {
  do.call(rbind, lapply(phi_grid, function(phi) {
    r <- compute_K(mu, phi)
    data.frame(mu = mu, phi = phi, K = r$K, iterations = r$iterations)
  }))
}))
results$abs_err <- abs(results$K - 1)
 
results$mu_f <- factor(results$mu, levels = c(2, 5, 10),
                        labels = c("mu==2", "mu==5", "mu==10"))
 
pal <- c("mu==2" = "#1b6ca8", "mu==5" = "#c9622b", "mu==10" = "#3a923a")
 
p1 <- ggplot(results, aes(x = phi, y = K, color = mu_f)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey50", linewidth = 0.4) +
  scale_x_log10() +
  scale_color_manual(values = pal, labels = scales::parse_format()) +
  labs(x = expression(phi ~ "(dispersion)"), y = expression(K(mu, phi)),
       title = "(a) Exact normalising constant", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
        plot.title = element_text(size = 12),
        panel.grid.minor = element_blank())
 
p2 <- ggplot(results, aes(x = phi, y = abs_err, color = mu_f)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey50", linewidth = 0.4) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = pal, labels = scales::parse_format()) +
  labs(x = expression(phi ~ "(dispersion)"), y = expression("|" * K(mu, phi) - 1 * "|"),
       title = "(b) Efron approximation error", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        plot.title = element_text(size = 12),
        panel.grid.minor = element_blank())

g <- arrangeGrob(p1, p2, ncol = 2)
ggsave("figures/dp_efron_error_R.png", g, width = 9.5, height = 4, dpi = 200)
cat("saved\n")

