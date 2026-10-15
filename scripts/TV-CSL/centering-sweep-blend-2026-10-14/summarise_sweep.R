# Summarise the centering sweep over every effect scale m that was run.
# Run from project root: Rscript scripts/TV-CSL/centering-sweep-blend-2026-10-14/summarise_sweep.R
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
D <- "scripts/TV-CSL/centering-sweep-blend-2026-10-14/res"
r <- bind_rows(lapply(list.files(D, "^rep_.*\\.csv$", full.names = TRUE), read.csv))
lab <- c("cox-linear-censored-only" = "marginal", "cox-time-varying-prop" = "time-varying",
         "cox-time-varying-oracle" = "oracle")
r$centering <- factor(lab[r$prop], levels = lab)

covr <- function(b, se, t) { ok <- is.finite(b) & is.finite(se); mean(abs(b[ok] - t[ok]) <= 1.96 * se[ok]) }
mse <- r |> distinct(n, m, centering, i, MSE, fit_secs) |> group_by(n, m, centering) |>
  summarise(reps = n(), MCSE = sd(MSE) / sqrt(n()), MSE = mean(MSE), fit_min = mean(fit_secs) / 60, .groups = "drop")
slopes <- r |> filter(coef_idx %in% 2:4) |> group_by(n, m, centering) |>
  summarise(bias = mean(beta - true), emp_sd = sd(beta), mean_se = mean(se_naive),
            cov_naive = covr(beta, se_naive, true), cov_sand = covr(beta, se_sandwich, true), .groups = "drop")
tab <- left_join(mse, slopes, by = c("n", "m", "centering")) |> arrange(n, m, centering)
cat("Slopes = beta_2..beta_4 (true value m). emp_sd / mean_se pooled over the three slopes.\n")
print(as.data.frame(tab |> select(n, m, centering, reps, MSE, MCSE, bias, emp_sd, mean_se, cov_naive, cov_sand, fit_min)),
      digits = 3, row.names = FALSE)

# Paired MSE differences on the same datasets (negative = time-varying / oracle better than marginal)
w <- r |> distinct(n, m, centering, i, MSE) |> pivot_wider(names_from = centering, values_from = MSE)
pd <- function(col) w |> filter(is.finite(.data[[col]]), is.finite(marginal)) |>
  group_by(n, m) |> summarise(reps = n(), diff = mean(.data[[col]] - marginal),
                              se = sd(.data[[col]] - marginal) / sqrt(n()), .groups = "drop") |>
  mutate(z = diff / se, contrast = paste(col, "- marginal"))
cat("\nPaired MSE difference vs marginal centering:\n")
print(as.data.frame(bind_rows(pd("time-varying"), if ("oracle" %in% names(w)) pd("oracle"))), digits = 3, row.names = FALSE)
