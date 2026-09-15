# plot-results.R
# Generates two summary plots from simulation results.
# Run from project root: Rscript scripts/TV-CSL/plot-results.R
#
# Output:
#   writeups/figures/plot-mse.pdf
#   writeups/figures/plot-coverage.pdf
#
# Data source: scripts/TV-CSL/results/ (via generate-plot-data.R → simulation-summary.R)

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr); library(stringr)
})

# Loads mse_long and cov_long from the consolidated results directory.
source("scripts/TV-CSL/generate-plot-data.R")

# ---- aesthetics ----
DGP_LABS <- c(
  "linear"     = "Linear η0  (linear spec.)",
  "non-linear" = "Non-linear η0  (complex spec.)"
)
METHOD_ORDER <- c(
  "S-Cox",
  "TV-CSL: marg. prop.",
  "TV-CSL: time-var. prop.",
  "TV-CSL: intercept-only"
)
METHOD_COLORS <- c(
  "S-Cox"                   = "#666666",
  "TV-CSL: marg. prop."     = "#0072B2",
  "TV-CSL: time-var. prop." = "#009E73",
  "TV-CSL: intercept-only"  = "#D55E00"
)
METHOD_SHAPES <- c(
  "S-Cox"                   = 17L,
  "TV-CSL: marg. prop."     = 16L,
  "TV-CSL: time-var. prop." = 15L,
  "TV-CSL: intercept-only"  = 18L
)
METHOD_LTYS <- c(
  "S-Cox"                   = "longdash",
  "TV-CSL: marg. prop."     = "solid",
  "TV-CSL: time-var. prop." = "solid",
  "TV-CSL: intercept-only"  = "dotdash"
)

base_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position   = "bottom",
    legend.title      = element_blank(),
    legend.key.width  = unit(1.8, "cm"),
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey92", color = NA),
    strip.text        = element_text(size = 10.5),
    plot.margin       = margin(4, 8, 4, 4)
  )

n_breaks <- c(200, 500, 1000, 2000)
n_labels <- c("200", "500", "1K", "2K")

factor_data <- function(d) {
  d |>
    mutate(
      Method = factor(Method, levels = METHOD_ORDER),
      DGP    = factor(DGP_LABS[DGP], levels = DGP_LABS)
    )
}

# ---- Plot 1: MSE ----
p_mse <- factor_data(mse_long) |>
  ggplot(aes(x = n, y = MSE, color = Method, shape = Method, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  scale_x_log10(breaks = n_breaks, labels = n_labels) +
  scale_y_log10(
    breaks = c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2),
    labels = c("0.02", "0.05", "0.1", "0.2", "0.5", "1", "2")
  ) +
  scale_color_manual(values = METHOD_COLORS, breaks = METHOD_ORDER) +
  scale_shape_manual(values = METHOD_SHAPES, breaks = METHOD_ORDER) +
  scale_linetype_manual(values = METHOD_LTYS, breaks = METHOD_ORDER) +
  facet_wrap(~DGP, ncol = 2, scales = "free_y") +
  labs(x = "Sample size  n", y = "MSE  (log scale)") +
  guides(color    = guide_legend(nrow = 2),
         shape    = guide_legend(nrow = 2),
         linetype = guide_legend(nrow = 2)) +
  base_theme

# ---- Plot 2: Coverage ----
p_cov <- factor_data(cov_long) |>
  ggplot(aes(x = n, y = cov, color = Method, shape = Method, linetype = Method)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", linewidth = 0.45) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  scale_x_log10(breaks = n_breaks, labels = n_labels) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    labels = paste0(seq(0, 100, 10), "%")
  ) +
  scale_color_manual(values = METHOD_COLORS, breaks = METHOD_ORDER) +
  scale_shape_manual(values = METHOD_SHAPES, breaks = METHOD_ORDER) +
  scale_linetype_manual(values = METHOD_LTYS, breaks = METHOD_ORDER) +
  facet_wrap(~DGP, ncol = 2) +
  labs(x = "Sample size  n",
       y = "Coverage of 95% CI  (naive SE)") +
  guides(color    = guide_legend(nrow = 2),
         shape    = guide_legend(nrow = 2),
         linetype = guide_legend(nrow = 2)) +
  base_theme +
  theme(panel.grid.major.y = element_line(color = "grey85"))

# ---- save ----
dir.create("writeups/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("writeups/figures/plot-mse.pdf",      p_mse, width = 6.5, height = 3.6,
       device = cairo_pdf)
ggsave("writeups/figures/plot-coverage.pdf", p_cov, width = 6.5, height = 3.6,
       device = cairo_pdf)
cat("Saved: writeups/figures/plot-mse.pdf\n")
cat("Saved: writeups/figures/plot-coverage.pdf\n")
