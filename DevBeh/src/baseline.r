#!/usr/bin/env Rscript
cat("Script started\n")
args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- dirname(normalizePath(script_path))
source(file.path(script_dir, "plots.r"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript script.R <input>")
}
input <- args[1]
cat("Input:", input, "\n")
run_BaselineH_plot(paste0(input,"Baseline_high/"))
run_BaselineL_plot(paste0(input,"Baseline_low/"))