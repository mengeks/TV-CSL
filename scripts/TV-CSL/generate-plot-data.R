# generate-plot-data.R
# Programmatic data extraction for plot-results.R.
# Sources simulation-summary.R (which reads scripts/TV-CSL/results/) and selects
# the right method × spec × DGP rows for the two main figures.
#
# Side effects: creates mse_long and cov_long in the calling environment.
#
# short_spec() in simulation-summary.R produces Spec values like:
#   TV_CSL:  "linear / linear-censored-only", "linear / intercept-only",
#            "linear / risk-set-adjusted", "linear / time-varying-prop",
#            "linear / time-varying-oracle",
#            "complex / linear-censored-only", etc.
#   S-Cox:   "linear", "complex"
#   Lasso:   "linear", "complex"
#
# Figure shows (per DGP):
#   S-Cox:               linear DGP → linear spec; non-linear DGP → complex spec
#   TV-CSL: marg. prop.  linear / linear-censored-only  (best η₀ spec per DGP)
#   TV-CSL: time-var.    linear / time-varying-prop      (best η₀ spec per DGP)
#   TV-CSL: intercept-only  linear / intercept-only      (best η₀ spec per DGP)
#
# "Best η₀ spec" convention: linear spec for linear DGP, complex spec for non-linear DGP.

suppressPackageStartupMessages(source("scripts/TV-CSL/simulation-summary.R"))
# agg and cov_agg are now available

# Display-name mapping function: (Method, Spec, eta_type) → plot label or NA to exclude
plot_label <- function(method, spec, eta_type) {
  if (method == "S-Cox") {
    if (eta_type == "linear"     & spec == "linear")  return("S-Cox")
    if (eta_type == "non-linear" & spec == "complex") return("S-Cox")
    return(NA_character_)
  }
  if (method == "TV_CSL") {
    eta_spec <- if (eta_type == "linear") "linear" else "complex"
    if (spec == paste0(eta_spec, " / linear-censored-only")) return("TV-CSL: marg. prop.")
    if (spec == paste0(eta_spec, " / time-varying-prop"))    return("TV-CSL: time-var. prop.")
    if (spec == paste0(eta_spec, " / intercept-only"))       return("TV-CSL: intercept-only")
    return(NA_character_)
  }
  NA_character_
}

# Apply labeling and filter
agg_labelled <- agg |>
  dplyr::rowwise() |>
  dplyr::mutate(PlotLabel = plot_label(Method, Spec, eta_type)) |>
  dplyr::ungroup() |>
  dplyr::filter(!is.na(PlotLabel))

cov_agg_labelled <- cov_agg |>
  dplyr::filter(coef_idx == 2) |>   # X.1 coefficient (first non-intercept; true β=1)
  dplyr::rowwise() |>
  dplyr::mutate(PlotLabel = plot_label(Method, Spec, eta_type)) |>
  dplyr::ungroup() |>
  dplyr::filter(!is.na(PlotLabel))

mse_long <- agg_labelled |>
  dplyr::select(DGP = eta_type, Method = PlotLabel, n, MSE)

cov_long <- cov_agg_labelled |>
  dplyr::select(DGP = eta_type, Method = PlotLabel, n, cov = cov_naive)
