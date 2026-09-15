# TV-CSL Simulation Study: Replication Guide

## Quick start (full replication from scratch)

```bash
bash scripts/TV-CSL/run-all.sh
```

This submits two SLURM array jobs (iterations 1–1000):

| Job name | Params file | Methods |
|---|---|---|
| `TV_CSL_main` | `params-main-methods.json` | S-Cox, S-Lasso, TV-CSL × {marg-prop, intercept-only, risk-set-adj., time-var-prop} |
| `TV_CSL_oracle` | `params-oracle-only.json` | TV-CSL × {time-var-oracle} (slower: numerical integration) |

Both jobs write to `scripts/TV-CSL/results/`. Data is auto-generated under `data/` if missing.

To submit a single custom run:
```bash
sbatch --job-name=MY_RUN scripts/TV-CSL/run-cluster.sh scripts/TV-CSL/params-main-methods.json
```

## Generate tables (after simulations complete)

```r
# From project root:
Rscript scripts/TV-CSL/generate-tables.R
```

Outputs LaTeX for all 6 tables to stdout:
- `tab:sim-main-mse`, `tab:sim-main-coverage` — condensed summary tables
- `tab:sim-linear-full`, `tab:sim-nonlinear-full` — full MSE tables by DGP
- `tab:cov-linear-full`, `tab:cov-nonlinear-full` — full coverage tables by DGP

## Generate figures

```r
Rscript scripts/TV-CSL/plot-results.R
```

Output: `writeups/figures/plot-mse.pdf` and `writeups/figures/plot-coverage.pdf`

## Check current results (without running full simulation)

```r
Rscript scripts/TV-CSL/simulation-summary.R
```

Prints MSE and 95% CI coverage aggregated from `scripts/TV-CSL/results/`.

## Compile writeup

```bash
cd writeups && pdflatex risk-set-adjust-tv-csl.tex
```

## Run unit tests

```bash
Rscript tests/test-tvcsl.R
```

---

## File layout

```
scripts/TV-CSL/
  run-cluster.sh              # Generic SLURM runner (takes params.json as $1)
  run-all.sh                  # Submits both SLURM jobs to replicate everything
  params-main-methods.json    # Config: S-Cox, S-Lasso, TV-CSL × 4 prop specs
  params-oracle-only.json     # Config: TV-CSL × oracle only (writes to same results/)
  simulation-summary.R        # Aggregates results/, prints MSE + coverage tables
  generate-tables.R           # Produces LaTeX for all 6 tables
  generate-plot-data.R        # Programmatic data extraction for figures
  plot-results.R              # Produces writeups/figures/*.pdf
  TV-CSL-runner.R             # Orchestration: reads JSON, runs methods, saves CSV
  time-varying-estimate.R     # Bridges TV-CSL-runner.R → R/old/TV-CSL.R
  results/                    # All simulation results (single source of truth)
  results-archive-2026-08/    # Archived: original main-methods results
  results-archive-tvprop-2026-08/   # Archived: time-varying-prop runs + bad oracle
  results-archive-oracle-2026-08/   # Archived: analytical oracle runs (partial)
```

## Variable role specification (JSON params)

TV-CSL-runner.R reads optional variable-role fields from `methods.TV_CSL` in the params JSON.
When omitted (or `null`), all `X.*` columns are used for all roles (default behavior):

```json
"TV_CSL": {
  "enabled": true,
  "prop_score_specs": ["cox-linear-censored-only"],
  "outcome_vars":     null,   // columns for η₀ (baseline hazard) model
  "treatment_vars":   null,   // columns for propensity score model
  "effect_modifiers": null    // columns for τ(x) HTE model
}
```

To use a subset, specify column names:
```json
"outcome_vars":     ["X.1", "X.2", "X.3"],
"treatment_vars":   ["X.2", "X.3"],
"effect_modifiers": ["X.1", "X.2", "X.3"]
```

## DGP parameters (current simulation)

| Parameter | Value |
|---|---|
| Covariates | p=3, X ~ Normal(0,I) |
| Sample sizes | n ∈ {200, 500, 1000, 2000} |
| Replications | R=1000 per cell |
| Cross-fitting folds | K=5 |
| Linear η₀ | η₀(x) = 2.5(x₁ + x₂/2 + x₃/3) |
| Non-linear η₀ | η₀(x) = −(3/4)σ(x₁)σ(x₂) |
| True HTE | τ₀(x) = x₁ + x₂ + x₃ |
| Censoring | λ_C = 0.1 (moderate) |
