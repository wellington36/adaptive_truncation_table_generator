suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(kableExtra)
})

# ---------------- plots ----------------
lev <- c("fixed_50", "fixed_200", "fixed_1000",
         "adaptive_1e-04", "adaptive_1e-07", "adaptive_1e-10")
lev <- lev[lev %in% unique(res$method)]
eps_of <- function(x) ifelse(grepl("adaptive", x),
                             as.numeric(sub("adaptive_", "", x)), NA_real_)

tab_main <- res %>%
  mutate(method = factor(method, lev)) %>%
  group_by(method) %>%
  summarise(
    `Coverage`          = mean(contains_ref, na.rm = TRUE),
    `p95 |error|`       = quantile(abs(ll_eval_err), .95, na.rm = TRUE),
    `Max |error|`       = max(abs(ll_eval_err), na.rm = TRUE),
    `% reps error>1e-6` = 100 * mean(abs(ll_eval_err) > 1e-6, na.rm = TRUE),
    `Max regret`        = max(regret, na.rm = TRUE),
    `Converged (%)`     = 100 * mean(conv == 0),
    .groups = "drop") %>%
  mutate(`Bound eps/2` = eps_of(as.character(method)) / 2, .after = 1)

kable(tab_main, format = "latex", booktabs = TRUE, digits = 3,
      format.args = list(scientific = TRUE),
      caption = "Evaluation error and regret, pooled over the 27 regimes.")

pal <- c(fixed_50 = "#d73027", fixed_200 = "#fc8d59", fixed_1000 = "#fdb863",
         `adaptive_1e-04` = "#9ecae1", `adaptive_1e-07` = "#4292c6",
         `adaptive_1e-10` = "#08306b")
line_pal <- c(fixed = "#e6784f", adaptive = "#2b6cb0")
FLOOR <- 1e-13

sm <- summarise_experiment1(res) %>%
  mutate(method = factor(method, lev),
         family = ifelse(startsWith(as.character(method), "adaptive"),
                         "adaptive", "fixed"),
         err    = pmax(p95_abs_ll_err, FLOOR),
         p_lab  = paste0("p_s = ", p_sent),
         R0_lab = paste0("R0 = ", R0)) %>%
  arrange(k, R0, p_sent, family, median_terms_total)

for (kk in c(.1, 1, 10)) {
  d <- sm %>% filter(k == kk)

  p <- ggplot(d, aes(median_terms_total, err)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = FLOOR * 1.5,
             fill = "grey85", alpha = .6) +
    geom_line(aes(group = family, colour = family),
              linewidth = 1.1, alpha = .55, show.legend = FALSE) +
    geom_point(aes(fill = method, shape = family),
               size = 4.5, colour = "white", stroke = .8) +
    facet_grid(p_lab ~ R0_lab) +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = c(FLOOR, NA)
    ) +
    scale_colour_manual(values = line_pal) +
    scale_fill_manual(values = pal, breaks = lev, name = NULL) +
    scale_shape_manual(values = c(fixed = 24, adaptive = 21), guide = "none") +
    guides(fill = guide_legend(nrow = 1, override.aes = list(shape = 21, size = 4))) +
    labs(x = "median total terms during optimisation",
         y = "p95 |log-lik error|",
         title = paste0("k = ", kk)) +
    theme_minimal(base_size = 21) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 15),  # Smaller method labels
      legend.key.width = unit(0.1, "cm"),  # Less space for each legend key
      legend.spacing.x = unit(0.2, "cm"),     # Less horizontal distance between items
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.1, "cm"),
      strip.text = element_text(face = "bold")
    )

  ggsave(sprintf("exp1_k%s.png", kk), p, width = 11, height = 7, dpi = 300)
}

grid27 <- expand_grid(R0 = c(.5, .7, .9), k = c(.1, 1, 10), p_sent = c(.01, .1, .5))
M <- 500; EPS <- 1e-10

iter_tab <- grid27 %>%
  rowwise() %>%
  mutate(
    n = tryCatch(adaptive_rho(R0, k, p_sent, EPS, M, max_n = 5e6L)$n,
                 error = function(e) NA_integer_),
    # log-lik error of a fixed cut-off J = 1000 relative to a 1e-12 reference
    err_J1000 = tryCatch({
      ref <- adaptive_rho(R0, k, p_sent, 1e-12, M, max_n = 5e6L)$log1m_mid
      abs(M * (log1m_rho_fixed(R0, k, p_sent, 1000L) - ref))
    }, error = function(e) NA_real_)
  ) %>% ungroup()

print(iter_tab, n = 27)

# ---- 9 x 3 LaTeX table: rows R0 x k, columns p_sent; bold if n > 1000 ----
wide <- iter_tab %>%
  mutate(cell = ifelse(n > 1000,
                       sprintf("\\textbf{%s}", format(n, big.mark = "{,}")),
                       format(n, big.mark = "{,}")),
         p = paste0("$p=", p_sent, "$")) %>%
  select(R0, k, p, cell) %>%
  pivot_wider(names_from = p, values_from = cell) %>%
  mutate(R0 = sprintf("$R_0=%s$", R0), k = sprintf("$k=%s$", k))
names(wide)[1:2] <- c("", "")

kable(wide, format = "latex", booktabs = TRUE, escape = FALSE, align = "llrrr",
      caption = "Terms needed by the adaptive rule ($\\varepsilon=10^{-10}$, $m=500$). Bold: more than 1000 terms.")

p_it <- iter_tab %>%
  ggplot(aes(factor(R0), factor(p_sent),
             fill = log10(n))) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(
    aes(
      label = format(n, big.mark = ","),
      fontface = ifelse(n > 1000, "bold", "plain"),
      colour = ifelse(n >= 1000, "white", "black")
    ),
    size = 5
  ) +
  facet_wrap(~ k, labeller = label_both) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    name = "Terms",
    breaks = log10(c(10, 100, 1000)),
    labels = c("10", "100", "1,000")
  ) +
  scale_colour_identity() +
  labs(
    x = expression(R[0]),
    y = expression(p[sent]),
    title = "Terms required for error ≤ 1e-10"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.spacing = unit(1.2, "lines"),
    strip.text = element_text(face = "bold", size = 14),
    axis.text = element_text(colour = "black")
  )

ggsave(
  "exp1_iterations.png",
  p_it,
  width = 12,
  height = 5,
  dpi = 300
)