# =============================================================
# 05_analyze_results.R
# Load MC results and produce summary tables + plots
# =============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

.args2 <- commandArgs(trailingOnly = FALSE)
.fa2   <- sub("--file=", "", .args2[grep("--file=", .args2)])
SCRIPT_DIR <- if (length(.fa2) && nchar(.fa2)) {
  normalizePath(dirname(.fa2), mustWork = FALSE)
} else {
  tryCatch(normalizePath(dirname(sys.frame(1)$ofile), mustWork = FALSE),
           error = function(e) getwd())
}
trailing_args <- commandArgs(trailingOnly = TRUE)
RESULTS_DIR <- if (length(trailing_args) >= 1 && nchar(trailing_args[1])) {
  normalizePath(trailing_args[1], mustWork = FALSE)
} else {
  file.path(SCRIPT_DIR, "results")
}

# ---- Load results ----
# Prefer individual files so partial runs are usable
iter_files <- sort(list.files(RESULTS_DIR, "mc_iter_.*\\.rds$", full.names = TRUE))

if (length(iter_files) == 0) {
  all_file <- file.path(RESULTS_DIR, "all_results.rds")
  if (!file.exists(all_file)) stop("No results found in ", RESULTS_DIR)
  all_res <- readRDS(all_file)
} else {
  all_res <- lapply(iter_files, readRDS)
}

message(sprintf("Loaded %d MC iterations", length(all_res)))

# ---- Extract metric data frame ----
METHODS <- c(
  csf          = "CSF",
  tlearner     = "T-learner (VT)",
  src1         = "S-learner (SRC1)",
  src2         = "SRC2",
  cox          = "Cox Interaction",
  ipcw_cf      = "IPCW-CF",
  oracle_lcsf  = "Oracle Linear CSF"
)
METRICS <- c("MSE", "RMSE", "Bias", "Cor", "SignError")

extract_df <- function(res_list, methods, metrics) {
  rows <- lapply(seq_along(res_list), function(k) {
    r <- res_list[[k]]
    if (inherits(r, "error") || is.null(r)) return(NULL)
    lapply(names(methods), function(m) {
      v <- r[[m]]
      if (is.null(v)) v <- setNames(rep(NA_real_, length(metrics)), metrics)
      data.frame(
        mc     = r$mc %||% k,
        method = m,
        as.list(v[metrics]),
        stringsAsFactors = FALSE
      )
    })
  })
  bind_rows(unlist(rows, recursive = FALSE))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

df_all <- extract_df(all_res, METHODS, METRICS)
df_all$method_label <- METHODS[df_all$method]
df_all$method_label <- factor(df_all$method_label, levels = unname(METHODS))

# ---- Summary table ----
tbl <- df_all %>%
  group_by(method_label) %>%
  summarise(
    n_ok         = sum(!is.na(MSE)),
    MSE_mean     = mean(MSE,       na.rm = TRUE),
    MSE_sd       = sd(MSE,         na.rm = TRUE),
    RMSE_mean    = mean(RMSE,      na.rm = TRUE),
    Bias_mean    = mean(Bias,      na.rm = TRUE),
    Cor_mean     = mean(Cor,       na.rm = TRUE),
    Cor_sd       = sd(Cor,         na.rm = TRUE),
    SignErr_mean = mean(SignError,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(MSE_mean)

cat("\n=== MSE Summary (lower is better) ===\n")
print(as.data.frame(tbl), digits = 4, row.names = FALSE)

# Save table
write.csv(tbl, file.path(RESULTS_DIR, "summary_table.csv"), row.names = FALSE)
message("Saved: ", file.path(RESULTS_DIR, "summary_table.csv"))

# ---- Plot 1: MSE boxplot ----
p1 <- ggplot(df_all, aes(x = method_label, y = MSE, fill = method_label)) +
  geom_boxplot(outlier.size = 0.8, show.legend = FALSE) +
  scale_x_discrete(guide = guide_axis(angle = 30)) +
  labs(title = "RMST-CATE MSE across 100 MC iterations",
       x = NULL, y = "MSE") +
  theme_bw(base_size = 12)

ggsave(file.path(RESULTS_DIR, "fig_mse_boxplot.pdf"), p1, width = 9, height = 5)
message("Saved: fig_mse_boxplot.pdf")

# ---- Plot 2: Correlation boxplot ----
p2 <- ggplot(df_all, aes(x = method_label, y = Cor, fill = method_label)) +
  geom_boxplot(outlier.size = 0.8, show.legend = FALSE) +
  scale_x_discrete(guide = guide_axis(angle = 30)) +
  labs(title = "Correlation(tau_hat, tau_true) across 100 MC iterations",
       x = NULL, y = "Pearson correlation") +
  theme_bw(base_size = 12)

ggsave(file.path(RESULTS_DIR, "fig_cor_boxplot.pdf"), p2, width = 9, height = 5)
message("Saved: fig_cor_boxplot.pdf")

# ---- Plot 3: tau_hat vs tau_true for first available MC iteration ----
first_res <- all_res[[which(sapply(all_res, function(r) !inherits(r, "error")))[1]]]
tau_true  <- first_res$tau_true

scatter_data <- bind_rows(lapply(names(METHODS), function(m) {
  th <- first_res[[paste0(m, "_hat")]]
  if (is.null(th)) return(NULL)
  data.frame(tau_true = tau_true, tau_hat = th,
             method = METHODS[m], stringsAsFactors = FALSE)
}))
scatter_data$method <- factor(scatter_data$method, levels = unname(METHODS))

p3 <- ggplot(scatter_data, aes(x = tau_true, y = tau_hat)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.6) +
  facet_wrap(~method, scales = "free_y") +
  labs(title = "Predicted vs True RMST-CATE (first MC iteration)",
       x = "True tau", y = "Estimated tau") +
  theme_bw(base_size = 10)

ggsave(file.path(RESULTS_DIR, "fig_scatter.pdf"), p3, width = 12, height = 8)
message("Saved: fig_scatter.pdf")

# ---- Oracle linear CSF beta coefficients ----
beta_list <- lapply(all_res, function(r) {
  if (inherits(r, "error") || is.null(r$oracle_lcsf_beta)) return(NULL)
  setNames(r$oracle_lcsf_beta, c("beta0", "beta1"))
})
beta_mat <- do.call(rbind, Filter(Negate(is.null), beta_list))
if (!is.null(beta_mat) && nrow(beta_mat) > 0) {
  cat("\n=== Oracle Linear CSF beta coefficients ===\n")
  cat(sprintf("  beta0 (intercept): mean=%.4f  sd=%.4f\n",
              mean(beta_mat[,1]), sd(beta_mat[,1])))
  cat(sprintf("  beta1 (X1 coef):   mean=%.4f  sd=%.4f\n",
              mean(beta_mat[,2]), sd(beta_mat[,2])))
}

message("Analysis complete.")
