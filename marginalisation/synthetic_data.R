source("experiments.R")

# Readable Experiment 1 plots: split by sentinel probability.
# Input: the raw tibble returned by run_simulation().
# This script is provided as code only and has not been run.

library(ggplot2)
library(dplyr)
library(purrr)

plot_experiment1_accuracy <- function(res,
                                      out_dir = ".",
                                      methods_to_show = c("adaptive", "fixed_20",
                                                          "fixed_50", "fixed_200",
                                                          "fixed_1000")) {
  required <- c("R0", "k", "p_sent", "replicate", "method", "abs_ll_error")
  missing <- setdiff(required, names(res))
  if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

  # Do not plot the reference as a competing method: its error against itself is zero.
  dat <- res %>%
    filter(method %in% methods_to_show, is.finite(abs_ll_error), abs_ll_error >= 0) %>%
    mutate(
      method = factor(method, levels = methods_to_show),
      R0 = factor(R0, levels = sort(unique(R0))),
      k = factor(k, levels = sort(unique(k))),
      p_sent = factor(p_sent, levels = sort(unique(p_sent)))
    )

  if (!nrow(dat)) stop("No rows remain after filtering to methods_to_show.")
  positives <- dat$abs_ll_error[dat$abs_ll_error > 0]
  if (!length(positives)) stop("All errors are zero; there is no positive log-scale range to plot.")

  # Log axes cannot display exact zero. Plot zeros at a disclosed floor below the
  # smallest positive error, and distinguish them with a separate point shape.
  zero_floor <- min(positives) / 3
  dat <- dat %>% mutate(
    error_for_plot = if_else(abs_ll_error == 0, zero_floor, abs_ll_error),
    is_zero = abs_ll_error == 0
  )
  y_limits <- c(zero_floor, max(dat$error_for_plot, na.rm = TRUE) * 1.15)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  sentinels <- levels(dat$p_sent)

  plots <- map(sentinels, function(ps) {
    d <- filter(dat, p_sent == ps)
    p <- ggplot(d, aes(x = method, y = error_for_plot, color = method)) +
      geom_boxplot(aes(group = method), width = 0.62, outlier.shape = NA,
                   alpha = 0.16, linewidth = 0.45) +
      geom_point(aes(shape = is_zero),
                 position = position_jitter(width = 0.10, height = 0),
                 size = 1.8, alpha = 0.65) +
      facet_grid(rows = vars(k), cols = vars(R0), labeller = label_both) +
      scale_y_log10(limits = y_limits) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 6),
                         labels = c(`FALSE` = "Positive error", `TRUE` = "Zero (shown at floor)")) +
      labs(
        title = paste0("Likelihood approximation accuracy (p_sent = ", ps, ")"),
        subtitle = paste0("Replicate errors versus the tighter adaptive reference; zero errors are placed at ",
                          format(zero_floor, scientific = TRUE), " for display only"),
        x = NULL,
        y = "Absolute log-likelihood error (log scale)",
        color = "Method",
        shape = "Observed error"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 25, hjust = 1),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold")
      )

    file_name <- file.path(out_dir, paste0("experiment1_accuracy_pSent_", ps, ".png"))
    ggsave(file_name, p, width = 9, height = 7, dpi = 320)
    p
  })
  names(plots) <- paste0("p_sent_", sentinels)
  plots
}

# Summarize only adaptive fits: average number of series terms over replicates
# for each parameter combination. CI describes Monte Carlo variation in the mean.
summarise_adaptive_terms <- function(res) {
  required <- c("R0", "k", "p_sent", "replicate", "method", "terms")
  missing <- setdiff(required, names(res))
  if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

  out <- res %>%
    filter(method == "adaptive", is.finite(terms)) %>%
    group_by(R0, k, p_sent) %>%
    summarise(
      n_replicates = n(),
      mean_terms = mean(terms),
      sd_terms = if (n() > 1L) sd(terms) else NA_real_,
      se_terms = if (n() > 1L) sd(terms) / sqrt(n()) else NA_real_,
      t_critical = if (n() > 1L) qt(0.975, df = n() - 1L) else NA_real_,
      ci_low = if (n() > 1L) max(0, mean(terms) - t_critical * se_terms) else NA_real_,
      ci_high = if (n() > 1L) mean(terms) + t_critical * se_terms else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      R0 = factor(R0, levels = sort(unique(R0))),
      k = factor(k, levels = sort(unique(k))),
      p_sent = factor(p_sent, levels = sort(unique(p_sent)))
    )
  out
}

# Create one mean adaptive-term plot for each sentinel probability.
plot_adaptive_terms <- function(term_summary, out_dir = ".") {
  required <- c("R0", "k", "p_sent", "mean_terms", "ci_low", "ci_high")
  missing <- setdiff(required, names(term_summary))
  if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  sentinels <- levels(factor(term_summary$p_sent))

  plots <- map(sentinels, function(ps) {
    d <- filter(term_summary, as.character(p_sent) == ps)
    p <- ggplot(d, aes(x = factor(R0), y = mean_terms, group = 1)) +
      geom_line(color = "#2878B5", linewidth = 0.55) +
      geom_point(color = "#176B87", size = 2.6) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                    width = 0.12, color = "#176B87", na.rm = TRUE) +
      facet_wrap(vars(k), nrow = 1, labeller = label_both) +
      scale_y_continuous(expand = expansion(mult = c(0.06, 0.12))) +
      labs(
        title = paste0("Adaptive series terms by model setting (p_sent = ", ps, ")"),
        subtitle = "Points show mean terms across replicates; bars show 95% t intervals for the mean",
        x = expression(R[0]),
        y = "Mean number of terms used by adaptive method"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold")
      )

    file_name <- file.path(out_dir, paste0("experiment1_adaptive_terms_pSent_", ps, ".png"))
    ggsave(file_name, p, width = 9, height = 4.8, dpi = 320)
    p
  })
  names(plots) <- paste0("p_sent_", sentinels)
  plots
}

# Example:
# Source the experiment code and run it with enough replicates to see distributions.
# results <- run_simulation(B = 100, m = 500, eps_ll = 1e-7,
#                           fixed_J = c(20L, 50L, 200L, 1000L))
# accuracy_plots <- plot_experiment1_accuracy(results, out_dir = "figures")
# accuracy_plots[[1]]  # display the first sentinel-probability plot
# terms <- summarise_adaptive_terms(results)
# write.csv(terms, "experiment1_adaptive_term_summary.csv", row.names = FALSE)
# term_plots <- plot_adaptive_terms(terms, out_dir = "figures")

# Interpretation:
# - Each output figure fixes p_sent, reducing the number of columns.
# - Rows facet by k; columns facet by R0.
# - Each box/point summarizes replicate-to-replicate log-likelihood error by method.
# - Lower is more accurate relative to the tighter adaptive reference.
# - The zero-error floor is only a display device; original values are unchanged.
# - Exclude "ref" from the method list because its error against itself is zero.
# - Keep a separate figure for runtime/terms: accuracy and computational cost are
#   distinct outcomes and should not be combined on one axis.
# - Adaptive term counts are averaged across replicates; report n_replicates too.
#   With one replicate, a confidence interval cannot be estimated (shown as NA).



sim <- run_simulation(
  B = 1L,
  m = 500L,
  eps_ll = 1e-12
)

accuracy_plots <- plot_experiment1_accuracy(sim, out_dir = "figures")
