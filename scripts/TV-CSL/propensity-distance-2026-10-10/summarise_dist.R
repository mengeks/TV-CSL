# Mean distance to the oracle a0, by m, evaluation set (events / risk sets) and propensity model.
# Run from project root: Rscript scripts/TV-CSL/propensity-distance-2026-10-10/summarise_dist.R
suppressPackageStartupMessages(library(dplyr))
r <- bind_rows(lapply(list.files("scripts/TV-CSL/propensity-distance-2026-10-10/res", "^dist_.*\\.csv$",
                                 full.names = TRUE), read.csv))
r$model <- factor(r$model, levels = c("marginal", "gam-additive", "gam-interaction"))
s <- r |> group_by(where, n, m, model) |>
  summarise(reps = n(), mae_se = sd(mae) / sqrt(n()), mae = mean(mae), rmse = mean(rmse), bias = mean(bias),
            .groups = "drop") |>
  group_by(where, n, m) |> mutate(ratio_to_marginal = mae / mae[model == "marginal"]) |> ungroup() |>
  arrange(where, n, m, model)
cat("MAE = mean |a_hat - a0|. ratio_to_marginal < 0.5 = candidate for a full TV-CSL run.\n")
print(as.data.frame(s |> select(where, n, m, model, reps, mae, mae_se, rmse, bias, ratio_to_marginal)),
      digits = 3, row.names = FALSE)
cat("\nMean secs per rep:", round(mean(r$secs), 1), "\n")
