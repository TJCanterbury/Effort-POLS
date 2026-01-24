#!/usr/bin/env Rscript
cat("Script started\n")
args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- dirname(normalizePath(script_path))

source(file.path(script_dir, "b_s.r"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript script.R <input>")
}
input <- args[1]
cat("Input:", input, "\n")

run_sigmacue_plot(paste0(input,"sigma_cue/"))
run_b_p_plot(paste0(input,"b_p/"))
run_c_q_plot(paste0(input,"c_q/"))
run_h_plot(paste0(input,"h/"))
run_theta_plot(paste0(input,"theta/"))
run_sigma_plot(paste0(input,"sigma/"))
run_cuecost_plot(paste0(input,"cue_cost/"))
run_b_s_plot(paste0(input,"b_s/"))