# True eta (this check) vs estimated eta (centering sweep), same datasets: m = 2, n = 4000.
# Run from project root: Rscript scripts/TV-CSL/eta-check-2026-10-10/summarise_eta.R
suppressPackageStartupMessages(library(dplyr))
tr <- bind_rows(lapply(list.files("scripts/TV-CSL/eta-check-2026-10-10/res", "^rep_", full.names = TRUE), read.csv))
reps <- sort(unique(tr$i))
es <- bind_rows(lapply(file.path("scripts/TV-CSL/centering-sweep-2026-10-10/res",
                                 sprintf("rep_n4000_m2_i%d.csv", reps)), read.csv)) |> mutate(eta = "estimated")
r <- bind_rows(tr, es) |> filter(prop %in% unique(tr$prop))
lab <- c("cox-linear-censored-only" = "marginal", "cox-time-varying-prop" = "time-varying (GAM)",
         "cox-time-varying-oracle" = "oracle")
r$centering <- lab[r$prop]
covr <- function(b, se, t) mean(abs(b - t) <= 1.96 * se, na.rm = TRUE)
out <- r |> group_by(centering, eta) |>
  summarise(reps = n_distinct(i), MSE = mean(MSE[coef_idx == 1]),
            slope_bias = mean((beta - true)[coef_idx %in% 2:4]),
            bias_se = sd((beta - true)[coef_idx %in% 2:4]) / sqrt(sum(coef_idx %in% 2:4)),
            coverage = covr(beta[coef_idx %in% 2:4], se_naive[coef_idx %in% 2:4], true[coef_idx %in% 2:4]),
            .groups = "drop") |> arrange(centering, eta)
print(as.data.frame(out), digits = 3, row.names = FALSE)
