# Pre-generate every dataset the main simulation reads (i = 1..1100: i trains, i + 100 tests),
# so the SLURM array never generates or reads a file another task is still writing.
# Uses the same ensure_data_exists() and JSON datagen params as TV-CSL-runner.R; existing files are kept.
# Run from project root: Rscript scripts/TV-CSL/pregenerate-data.R [params.json]
suppressPackageStartupMessages(source("scripts/TV-CSL/TV-CSL-runner.R"))
json_file <- commandArgs(TRUE)[1]
if (is.na(json_file)) json_file <- "scripts/TV-CSL/params-main-methods.json"
config <- fromJSON(json_file)
R <- ifelse(is.null(config$R), 1000, config$R)
t0 <- Sys.time()
for (eta_type in c("linear", "non-linear")) for (n in c(200, 500, 1000, 2000)) for (i in seq_len(R + 100)) {
  suppressMessages(invisible(capture.output(
    ensure_data_exists(i, n, eta_type, "linear", config$datagen))))
}
cat(sprintf("Done: %d datasets checked in %.1f min\n", 2 * 4 * (R + 100),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
