suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
library(patchwork)
})

loci_labels <- list(
  u_base = bquote("Responsiveness effort threshold "(u)),
  rho    = bquote("Baseline effort " (rho)),
  nu     = bquote("POLS effort bias "(nu)),
  gamma  = bquote("Social information effort bias "(gamma)),
  lambda = bquote("Responsiveness to partner effort "(lambda)),
  c      = bquote("Baseline observation rate "(c)),
  m      = bquote("POLS observation bias "(m))
)

run_trait_plot <- function(
  path,
  x_var,
  x_label,
  out_name,
  y_limits = c(-1, 1)
) {

  readfile <- read.csv(paste0(path, "summaries.csv"))

  # tidy evaluation
  x_sym <- rlang::sym(x_var)

  df_long <- readfile %>%
    dplyr::select(!!x_sym, u_base, rho, nu, gamma, lambda, c, m) %>%
    tidyr::gather(Loci, Value, -!!x_sym)

  x_limits <- range(readfile[[x_var]])

  p <- ggplot(df_long, aes(
      x = !!x_sym,
      y = Value,
      color = Loci,
      group = Loci,
      fill  = Loci
    )
  ) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
    theme_classic(base_size = 12) +
    xlab(x_label) +
    theme(
      panel.grid.major = element_line(colour = "grey80", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey90", linewidth = 0.2)
    ) +
    scale_color_discrete(name = "Loci", labels = loci_labels) +
    scale_fill_discrete(name = "Loci", labels = loci_labels) +
    coord_cartesian(xlim = x_limits, ylim = y_limits)

  pdf(paste0(path, "../", out_name, ".pdf"), width = 8, height = 6)
  print(p)
  dev.off()

  file.copy(
    paste0(path, "../", out_name, ".pdf"),
    paste0(path, "../../", out_name, ".pdf"),
    overwrite = TRUE
  )
}

run_b_f_plot <- function(path) {
  run_trait_plot(
    path,
    x_var   = "b_f",
    x_label = bquote("Independence of offspring "(b[f])),
    out_name = "b_f"
  )
}

run_h_plot <- function(path) {
  run_trait_plot(
    path,
    x_var   = "h",
    x_label = bquote("Hawkishness of fast pace-of-life individuals "(eta)),
    out_name = "h"
  )
}

run_theta_plot <- function(path) {
  run_trait_plot(
    path,
    x_var   = "theta",
    x_label = bquote("Slope of mortality for pace-of-life individuals "(theta)),
    out_name = "theta"
  )
}

run_sigma_plot <- function(path) {
  run_trait_plot(
    path,
    x_var   = "sigma0",
    x_label = bquote("Standard deviation of pace-of-life phenotype "(sigma)),
    out_name = "sigma0"
  )
}

run_sigmacue_plot <- function(path) {
  run_trait_plot(
    path,
    x_var   = "sigma_cue",
    x_label = bquote("Standard deviation of social cue"),
    out_name = "sigma_cue"
  )
}

run_cuecost_plot <- function(path) {
  run_trait_plot(
    path,
    x_var   = "c_v",
    x_label = bquote("Cost of social cue "(c[v])),
    out_name = "cue_cost"
  )
}

run_b_s_plot <- function(path) {
  run_trait_plot(
    path,
    x_var   = "b_s",
    x_label = bquote("Baseline survival rate "(b[s])),
    out_name = "b_s"
  )
}

run_b_p_plot <- function(path) {
  run_trait_plot(
    path,
    x_var   = "b_p",
    x_label = bquote("Brood predation risk "(b[p])),
    out_name = "b_p"
  )
}

run_divorce_plot <- function(path) {
  run_trait_plot(
    path,
    x_var   = "div_rate",
    x_label = bquote("Divorce rate"),
    out_name = "div_rate"
  )
}
