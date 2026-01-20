#!/usr/bin/env Rscript
cat("Script started\n")
args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- dirname(normalizePath(script_path))
library(akima)
source(file.path(script_dir, "b_s.r"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript script.R <input>")
}
input <- args[1]
cat("Input:", input, "\n")
run_Baseline_plot(paste0(input,"Baseline/"))