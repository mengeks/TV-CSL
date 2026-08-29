suppressPackageStartupMessages({
  library(vroom); library(dplyr); library(stringr)
})

# Edit this vector to select which results directories to include.
# Comment out any dirs you want to exclude; non-existent dirs are silently skipped.
results_dirs <- c(
  "scripts/TV-CSL/results",               # main 4 methods (marg-prop, intercept-only, etc.)
  "scripts/TV-CSL/results_new_props",     # old name for time-varying-prop/oracle runs
  "scripts/TV-CSL/results_time_varying_prop"  # new name (after cluster re-run)
)
results_dirs <- results_dirs[dir.exists(results_dirs)]
cat("Reading from dirs:\n"); cat(" ", results_dirs, sep = "\n  "); cat("\n")

# Normalize old spec names → current names
normalize_spec <- function(x) {
  x <- str_replace_all(x, "cox-direct-event",             "cox-time-varying-prop")
  x <- str_replace_all(x, "cox-risk-set-adjusted-oracle", "cox-time-varying-oracle")
  x
}

# ---- Result CSVs (batched loading) ----
f_res <- unlist(lapply(results_dirs, function(d)
  list.files(d, pattern = "result-iteration", recursive = TRUE, full.names = TRUE)))
f_res <- f_res[!grepl("_inference", f_res)]
cat("Result files:", length(f_res), "\n")

batch_size <- 200
batches <- split(f_res, ceiling(seq_along(f_res) / batch_size))
res_list <- lapply(batches, function(b) {
  d <- vroom(b, id = "path", show_col_types = FALSE, progress = FALSE)
  d$Specification <- normalize_spec(d$Specification)
  d$eta_type <- str_extract(d$path, "eta-([a-z-]+)_HTE") |> str_remove("eta-") |> str_remove("_HTE")
  d$n        <- as.integer(str_extract(d$path, "_n-(\\d+)/") |> str_remove_all("_n-|/"))
  d
})
res <- bind_rows(res_list)

short_spec <- function(method, spec) {
  s <- as.character(spec)
  ifelse(method == "TV_CSL",
    paste0(ifelse(grepl("regressor-spec-linear", s), "linear", "complex"),
           " / ", str_extract(s, "cox-[^\"]+") |> str_remove("cox-")),
    str_extract(s, "regressor-spec-[^_]+") |> str_remove("regressor-spec-"))
}
res$Spec <- short_spec(res$Method, res$Specification)

agg <- res |>
  group_by(Method, Spec, eta_type, n) |>
  summarise(MSE    = mean(MSE_Estimate, na.rm = TRUE),
            MCSE   = sd(MSE_Estimate,   na.rm = TRUE) / sqrt(n()),
            Time_s = mean(Time_Taken,   na.rm = TRUE),
            Runs   = n(), .groups = "drop") |>
  arrange(eta_type, n, Method, Spec)

cat("\n=== MSE Summary ===\n")
for (et in c("linear", "non-linear")) {
  cat(sprintf("\neta = %s\n", et))
  cat(sprintf("%-10s %-50s %5s  %7s (%6s)  %8s  [runs]\n",
              "Method", "Spec", "n", "MSE", "MCSE", "Time(s)"))
  cat(strrep("-", 95), "\n")
  sub <- filter(agg, eta_type == et)
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    cat(sprintf("%-10s %-50s %5d  %7.4f (%6.4f)  %8.1f  [%d]\n",
                r$Method, r$Spec, r$n, r$MSE, r$MCSE, r$Time_s, r$Runs))
  }
}

# ---- Inference CSVs ----
f_inf <- unlist(lapply(results_dirs, function(d)
  list.files(d, pattern = "_inference\\.csv$", recursive = TRUE, full.names = TRUE)))
cat("\n\nInference files:", length(f_inf), "\n")

batches_inf <- split(f_inf, ceiling(seq_along(f_inf) / batch_size))
inf_list <- lapply(batches_inf, function(b) {
  d <- vroom(b, id = "path", show_col_types = FALSE, progress = FALSE)
  d$Specification <- normalize_spec(d$Specification)
  d$eta_type <- str_extract(d$path, "eta-([a-z-]+)_HTE") |> str_remove("eta-") |> str_remove("_HTE")
  d$n        <- as.integer(str_extract(d$path, "_n-(\\d+)/") |> str_remove_all("_n-|/"))
  d
})
inf <- bind_rows(inf_list)

# True beta for HTE_type="linear", p=3: (intercept=0, X.1=1, X.2=1, X.3=1)
inf$true_beta <- c(0, 1, 1, 1)[inf$coef_idx]
inf$Spec <- short_spec(inf$Method, inf$Specification)

cov_fn <- function(beta, se, true_b) {
  ok <- is.finite(se) & is.finite(beta)
  if (!any(ok)) return(NA_real_)
  mean((beta[ok] - 1.96 * se[ok] <= true_b[ok]) & (true_b[ok] <= beta[ok] + 1.96 * se[ok]))
}

cov_agg <- inf |>
  group_by(Method, Spec, eta_type, n, coef_idx) |>
  summarise(
    cov_naive    = cov_fn(beta, se_naive,    true_beta),
    cov_sandwich = cov_fn(beta, se_sandwich, true_beta),
    Runs         = n(),
    .groups      = "drop"
  ) |>
  arrange(eta_type, n, Method, Spec, coef_idx)

cat("\n=== 95% CI Coverage ===\n")
fmt_cov <- function(x) ifelse(is.na(x), "  NA  ", sprintf("%.3f", x))

for (et in c("linear", "non-linear")) {
  cat(sprintf("\neta = %s\n", et))
  cat(sprintf("%-10s %-50s %5s  coef  %-7s  %-7s  [runs]\n",
              "Method", "Spec", "n", "cov_naive", "cov_sand"))
  cat(strrep("-", 95), "\n")
  sub    <- filter(cov_agg, eta_type == et)
  prev   <- ""
  for (i in seq_len(nrow(sub))) {
    r   <- sub[i, ]
    key <- paste(r$Method, r$Spec, r$n)
    if (key != prev) cat("\n")
    prev <- key
    cat(sprintf("%-10s %-50s %5d  [β%d]  %-7s  %-7s  [%d]\n",
                r$Method, r$Spec, r$n, r$coef_idx,
                fmt_cov(r$cov_naive), fmt_cov(r$cov_sandwich), r$Runs))
  }
}
