# generate-tables.R
# Generates all simulation tables for the writeup from scripts/TV-CSL/results/.
# Run from project root: Rscript scripts/TV-CSL/generate-tables.R
#
# Output: LaTeX table code printed to stdout for each of the 6 tables:
#   Full (4): tab:sim-linear-full, tab:sim-nonlinear-full,
#             tab:cov-linear-full,  tab:cov-nonlinear-full
#   Main (2): tab:sim-main-mse,    tab:sim-main-coverage
#
# The coverage tables average over coef_idx 2:4 (β₂, β₃, β₄; true value = 1).

suppressPackageStartupMessages(source("scripts/TV-CSL/simulation-summary.R"))
suppressPackageStartupMessages(library(tidyr))

# ---- Row-label helpers -------------------------------------------------------

prop_label <- function(spec) {
  dplyr::case_when(
    grepl("linear-censored-only", spec) ~ "marg.\\ prop.",
    grepl("intercept-only",       spec) ~ "intercept-only",
    grepl("risk-set-adjusted",    spec) ~ "risk-set-adj.",
    grepl("time-varying-prop",    spec) ~ "time-var.\\ prop.",
    grepl("time-varying-oracle",  spec) ~ "time-var.\\ oracle",
    TRUE ~ spec
  )
}

eta_label <- function(spec) {
  dplyr::if_else(grepl("^linear", spec), "linear", "complex")
}

# ---- Full MSE table ----------------------------------------------------------

make_full_mse_table <- function(et, label, caption) {
  sub <- agg |>
    dplyr::filter(eta_type == et) |>
    dplyr::mutate(
      MethodLabel = dplyr::case_when(
        Method == "TV_CSL" ~ "TV-CSL",
        Method == "S-Cox"  ~ "S-Cox",
        Method == "Lasso"  ~ "S-Lasso",
        TRUE ~ Method
      ),
      EtaSpec  = eta_label(Spec),
      PropSpec = prop_label(Spec)
    ) |>
    dplyr::select(MethodLabel, EtaSpec, PropSpec, n, MSE) |>
    tidyr::pivot_wider(names_from = n, values_from = MSE,
                       names_prefix = "n_") |>
    dplyr::arrange(MethodLabel, EtaSpec, PropSpec)

  fmt <- function(x) sprintf("%.3f", x)

  lines <- c(
    "\\begin{table}[h]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{llcccc}",
    "\\toprule",
    "Method & $\\eta_0$ / prop.\\ spec & $n=200$ & $n=500$ & $n=1000$ & $n=2000$ \\\\",
    "\\midrule"
  )

  prev_method <- ""
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    if (r$MethodLabel != prev_method && prev_method != "") {
      lines <- c(lines, "\\midrule")
    }
    prev_method <- r$MethodLabel

    if (r$MethodLabel %in% c("S-Cox", "S-Lasso")) {
      row_label <- paste0(r$MethodLabel, " & ", r$EtaSpec)
    } else {
      if (!is.na(r$PropSpec) && nchar(r$PropSpec) > 0) {
        row_label <- paste0(r$MethodLabel, " & ", r$EtaSpec, " / ", r$PropSpec)
      } else {
        row_label <- paste0(r$MethodLabel, " & ", r$EtaSpec)
      }
    }

    lines <- c(lines, paste0(
      row_label, " & ",
      fmt(r$n_200), " & ", fmt(r$n_500), " & ",
      fmt(r$n_1000), " & ", fmt(r$n_2000), " \\\\"
    ))
  }

  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  paste(lines, collapse = "\n")
}

# ---- Full coverage table -----------------------------------------------------

make_full_cov_table <- function(et, best_eta_spec, label, caption) {
  # best_eta_spec: "linear" for linear DGP, "complex" for non-linear DGP (TV-CSL)
  # S-Cox/Lasso: always use the best available spec

  cov_avg <- cov_agg |>
    dplyr::filter(coef_idx %in% 2:4) |>
    dplyr::filter(eta_type == et) |>
    dplyr::group_by(Method, Spec, eta_type, n) |>
    dplyr::summarise(
      cov_naive    = mean(cov_naive,    na.rm = TRUE),
      cov_sandwich = mean(cov_sandwich, na.rm = TRUE),
      Runs         = mean(Runs),
      .groups = "drop"
    )

  sub <- cov_avg |>
    dplyr::mutate(
      MethodLabel = dplyr::case_when(
        Method == "TV_CSL" ~ "TV-CSL",
        Method == "S-Cox"  ~ "S-Cox",
        Method == "Lasso"  ~ "S-Lasso",
        TRUE ~ Method
      ),
      EtaSpec  = eta_label(Spec),
      PropSpec = prop_label(Spec)
    ) |>
    dplyr::filter(
      (MethodLabel %in% c("S-Cox", "S-Lasso") & EtaSpec == "linear") |
      (MethodLabel == "TV-CSL" & EtaSpec == best_eta_spec)
    ) |>
    dplyr::select(MethodLabel, EtaSpec, PropSpec, n, cov_naive, cov_sandwich) |>
    tidyr::pivot_wider(
      names_from = n,
      values_from = c(cov_naive, cov_sandwich),
      names_sep = "_n"
    ) |>
    dplyr::arrange(MethodLabel, EtaSpec, PropSpec)

  fmt <- function(x) ifelse(is.na(x), "---  ", sprintf("%.3f", x))

  lines <- c(
    "\\begin{table}[h]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{llcccc}",
    "\\toprule",
    "Method & prop.\\ spec & $n=200$ & $n=500$ & $n=1000$ & $n=2000$ \\\\",
    "\\midrule"
  )

  prev_method <- ""; prev_se <- ""
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    if (r$MethodLabel != prev_method && prev_method != "") lines <- c(lines, "\\midrule")
    prev_method <- r$MethodLabel

    if (r$MethodLabel %in% c("S-Cox", "S-Lasso")) {
      lines <- c(lines, paste0(
        r$MethodLabel, " / S-Cox & --- & ",
        fmt(r$cov_naive_n200), " & ", fmt(r$cov_naive_n500), " & ",
        fmt(r$cov_naive_n1000), " & ", fmt(r$cov_naive_n2000), " \\\\"
      ))
    } else {
      lines <- c(lines, paste0(
        "TV-CSL (naive)    & ", r$PropSpec, " & ",
        fmt(r$cov_naive_n200), " & ", fmt(r$cov_naive_n500), " & ",
        fmt(r$cov_naive_n1000), " & ", fmt(r$cov_naive_n2000), " \\\\"
      ))
      lines <- c(lines, paste0(
        "TV-CSL (sandwich) & ", r$PropSpec, " & ",
        fmt(r$cov_sandwich_n200), " & ", fmt(r$cov_sandwich_n500), " & ",
        fmt(r$cov_sandwich_n1000), " & ", fmt(r$cov_sandwich_n2000), " \\\\"
      ))
    }
  }

  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  paste(lines, collapse = "\n")
}

# ---- Main MSE table (both DGPs side-by-side) ---------------------------------

make_main_mse_table <- function() {
  # S-Cox: linear spec for linear DGP, complex spec for non-linear DGP.
  # TV-CSL: linear spec for linear DGP, complex spec for non-linear DGP.
  # Show: S-Cox, TV-CSL × {marg-prop, intercept-only, time-var-prop, time-var-oracle}.

  select_main <- function(et) {
    eta_spec <- if (et == "linear") "linear" else "complex"
    bind_rows(
      agg |> dplyr::filter(Method == "S-Cox", eta_type == et, Spec == eta_spec),
      agg |> dplyr::filter(Method == "TV_CSL", eta_type == et,
                           Spec == paste0(eta_spec, " / linear-censored-only")),
      agg |> dplyr::filter(Method == "TV_CSL", eta_type == et,
                           Spec == paste0(eta_spec, " / intercept-only")),
      agg |> dplyr::filter(Method == "TV_CSL", eta_type == et,
                           Spec == paste0(eta_spec, " / time-varying-prop")),
      agg |> dplyr::filter(Method == "TV_CSL", eta_type == et,
                           Spec == paste0(eta_spec, " / time-varying-oracle"))
    ) |>
      dplyr::mutate(
        RowLabel = dplyr::case_when(
          Method == "S-Cox" ~ "S-Cox",
          grepl("linear-censored-only", Spec) ~ "TV-CSL: marg.\\ prop.",
          grepl("intercept-only",       Spec) ~ "TV-CSL: intercept-only",
          grepl("time-varying-prop",    Spec) ~ "TV-CSL: time-var.\\ prop.",
          grepl("time-varying-oracle",  Spec) ~ "TV-CSL: time-var.\\ oracle",
          TRUE ~ NA_character_
        )
      ) |>
      dplyr::filter(!is.na(RowLabel)) |>
      dplyr::select(RowLabel, n, MSE) |>
      tidyr::pivot_wider(names_from = n, values_from = MSE, names_prefix = "n_")
  }

  lin <- select_main("linear")
  nonlin <- select_main("non-linear")

  stopifnot(identical(lin$RowLabel, nonlin$RowLabel))

  fmt <- function(x) sprintf("%.3f", x)

  lines <- c(
    "\\begin{table}[h]",
    "\\centering",
    "\\caption{MSE summary. Linear $\\eta_0$: linear $\\eta_0$ spec; non-linear $\\eta_0$: complex $\\eta_0$ spec. Lower is better.}",
    "\\label{tab:sim-main-mse}",
    "\\begin{tabular}{lcccccccc}",
    "\\toprule",
    " & \\multicolumn{4}{c}{Linear $\\eta_0$} & \\multicolumn{4}{c}{Non-linear $\\eta_0$} \\\\",
    "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
    "Method & $n{=}200$ & $n{=}500$ & $n{=}1000$ & $n{=}2000$ & $n{=}200$ & $n{=}500$ & $n{=}1000$ & $n{=}2000$ \\\\",
    "\\midrule"
  )

  prev_group <- ""
  for (i in seq_len(nrow(lin))) {
    group <- if (grepl("S-Cox", lin$RowLabel[i])) "scox" else "tvcsl"
    if (group != prev_group && prev_group != "") lines <- c(lines, "\\midrule")
    prev_group <- group

    lines <- c(lines, paste0(
      lin$RowLabel[i], " & ",
      fmt(lin$n_200[i]),  " & ", fmt(lin$n_500[i]),  " & ",
      fmt(lin$n_1000[i]), " & ", fmt(lin$n_2000[i]), " & ",
      fmt(nonlin$n_200[i]),  " & ", fmt(nonlin$n_500[i]),  " & ",
      fmt(nonlin$n_1000[i]), " & ", fmt(nonlin$n_2000[i]), " \\\\"
    ))
  }

  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  paste(lines, collapse = "\n")
}

# ---- Main coverage table (both DGPs side-by-side, naive SE only) ------------

make_main_cov_table <- function() {
  select_main_cov <- function(et) {
    eta_spec <- if (et == "linear") "linear" else "complex"
    cov_agg |>
      dplyr::filter(coef_idx %in% 2:4, eta_type == et) |>
      dplyr::group_by(Method, Spec, n) |>
      dplyr::summarise(cov_naive = mean(cov_naive, na.rm = TRUE), .groups = "drop") |>
      dplyr::filter(
        (Method == "S-Cox"  & Spec == "linear") |
        (Method == "TV_CSL" & Spec == paste0(eta_spec, " / linear-censored-only")) |
        (Method == "TV_CSL" & Spec == paste0(eta_spec, " / intercept-only")) |
        (Method == "TV_CSL" & Spec == paste0(eta_spec, " / time-varying-prop")) |
        (Method == "TV_CSL" & Spec == paste0(eta_spec, " / time-varying-oracle"))
      ) |>
      dplyr::mutate(
        RowLabel = dplyr::case_when(
          Method == "S-Cox" ~ "S-Cox",
          grepl("linear-censored-only", Spec) ~ "TV-CSL: marg.\\ prop.",
          grepl("intercept-only",       Spec) ~ "TV-CSL: intercept-only",
          grepl("time-varying-prop",    Spec) ~ "TV-CSL: time-var.\\ prop.",
          grepl("time-varying-oracle",  Spec) ~ "TV-CSL: time-var.\\ oracle",
          TRUE ~ NA_character_
        )
      ) |>
      dplyr::filter(!is.na(RowLabel)) |>
      dplyr::select(RowLabel, n, cov_naive) |>
      tidyr::pivot_wider(names_from = n, values_from = cov_naive, names_prefix = "n_")
  }

  lin <- select_main_cov("linear")
  nonlin <- select_main_cov("non-linear")

  stopifnot(identical(lin$RowLabel, nonlin$RowLabel))

  fmt <- function(x) sprintf("%.3f", x)

  lines <- c(
    "\\begin{table}[h]",
    "\\centering",
    "\\caption{Empirical 95\\% CI coverage (naive SE), averaged over $\\beta_2,\\beta_3,\\beta_4$. Linear $\\eta_0$: linear spec; non-linear $\\eta_0$: complex spec.}",
    "\\label{tab:sim-main-coverage}",
    "\\begin{tabular}{lcccccccc}",
    "\\toprule",
    " & \\multicolumn{4}{c}{Linear $\\eta_0$} & \\multicolumn{4}{c}{Non-linear $\\eta_0$} \\\\",
    "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
    "Method & $n{=}200$ & $n{=}500$ & $n{=}1000$ & $n{=}2000$ & $n{=}200$ & $n{=}500$ & $n{=}1000$ & $n{=}2000$ \\\\",
    "\\midrule"
  )

  prev_group <- ""
  for (i in seq_len(nrow(lin))) {
    group <- if (grepl("S-Cox", lin$RowLabel[i])) "scox" else "tvcsl"
    if (group != prev_group && prev_group != "") lines <- c(lines, "\\midrule")
    prev_group <- group

    lines <- c(lines, paste0(
      lin$RowLabel[i], " & ",
      fmt(lin$n_200[i]),  " & ", fmt(lin$n_500[i]),  " & ",
      fmt(lin$n_1000[i]), " & ", fmt(lin$n_2000[i]), " & ",
      fmt(nonlin$n_200[i]),  " & ", fmt(nonlin$n_500[i]),  " & ",
      fmt(nonlin$n_1000[i]), " & ", fmt(nonlin$n_2000[i]), " \\\\"
    ))
  }

  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  paste(lines, collapse = "\n")
}

# ---- Output all tables -------------------------------------------------------

cat("\n\n% ===== tab:sim-linear-full =====\n\n")
cat(make_full_mse_table(
  "linear",
  "tab:sim-linear-full",
  "MSE under linear $\\eta_0$. All $\\eta_0$ and propensity specifications. Lower is better."
))

cat("\n\n% ===== tab:sim-nonlinear-full =====\n\n")
cat(make_full_mse_table(
  "non-linear",
  "tab:sim-nonlinear-full",
  "MSE under non-linear $\\eta_0$. All $\\eta_0$ and propensity specifications."
))

cat("\n\n% ===== tab:cov-linear-full =====\n\n")
cat(make_full_cov_table(
  "linear", "linear",
  "tab:cov-linear-full",
  "Empirical 95\\% CI coverage under linear $\\eta_0$ (linear $\\eta_0$ spec). Averaged over $\\beta_2,\\beta_3,\\beta_4$."
))

cat("\n\n% ===== tab:cov-nonlinear-full =====\n\n")
cat(make_full_cov_table(
  "non-linear", "complex",
  "tab:cov-nonlinear-full",
  "Empirical 95\\% CI coverage under non-linear $\\eta_0$ (complex $\\eta_0$ spec). Averaged over $\\beta_2,\\beta_3,\\beta_4$."
))

cat("\n\n% ===== tab:sim-main-mse =====\n\n")
cat(make_main_mse_table())

cat("\n\n% ===== tab:sim-main-coverage =====\n\n")
cat(make_main_cov_table())
