suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
library(patchwork)
library(akima)
library(fields)
library(ggpubr)
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

trait_labels <- list(
  mean_faster_effort = bquote("Effort of fast individuals"),
  mean_slower_effort    = bquote("Effort of slow individuals"),
  effort_dif    = bquote("Asymmetries in effort within pairs"),
  mean_f_s_d    = bquote("Effort distribution in fast-slow pairs"),
  mean_f_s_w    = bquote("fast-slow pairs"),
  f_w    = bquote("fast-only pairs"),
  s_w    = bquote("slow-only pairs"),
  mean_fast_h     = bquote("Uncertainty of fast individuals"),
  mean_slow_h  = bquote("Uncertainty of slow individuals")
)

run_BaselineL_plot <- function(path) {
  
  readfile <- read.csv(paste0(path, "summaries.csv"))

  df <- readfile %>%
    mutate(r=u1/(u1+u2)) 
  
  
  # Create grid
  gx <- seq(-1, 1, length = 20)
  gy <- seq(-1, 1, length = 20)


  # Predation
  pdf(paste0(path, "../", "baselineL_pred", ".pdf"), width = 8, height = 8)
  # Interpolate using akima-style linear interpolation
  interp <- akima::interp(df$qm, df$qf, 1-df$pred, xo = gx, yo = gy, duplicate = "mean")
  # Heatmap
  image.plot(
    interp$x, interp$y, interp$z,
    col  = hcl.colors(100, "Inferno"),
    zlim = c(0,1),
    xlab = expression(q[m]),
    ylab = expression(q[f]),
    legend.lab = expression("Nest defence")
  )
  # Data support
  # points(df$qm, df$qf, pch = 16, cex = 0.2, col = rgb(1,1,1,0.3))
  dev.off()

  # Effort
  pdf(paste0(path, "../", "baselineL_u", ".pdf"), width = 8, height = 8)
  # Interpolate using akima-style linear interpolation
  interp <- akima::interp(df$qm, df$qf, df$r, xo = gx, yo = gy, duplicate = "mean")
  # Heatmap
  image.plot(
    interp$x, interp$y, interp$z,
    col  = hcl.colors(100, "Inferno"),
    zlim = c(0,1),
    xlab = expression(q[m]),
    ylab = expression(q[f]),
    legend.lab = expression("male:female effort ratio"),
    cex=1.5
  )
  # Data support
  # points(df$qm, df$qf, pch = 16, cex = 0.2, col = rgb(1,1,1,0.3))
  dev.off()

  file.copy(
    paste0(path, "../", "baselineL_u", ".pdf"),
    paste0(path, "../../", "baselineL_u", ".pdf"),
    overwrite = TRUE
  )
  file.copy(
    paste0(path, "../", "baselineL_pred", ".pdf"),
    paste0(path, "../../", "baselineL_pred", ".pdf"),
    overwrite = TRUE
  )
}

run_BaselineH_plot <- function(path) {
  
  readfile <- read.csv(paste0(path, "summaries.csv"))

  df <- readfile %>%
    mutate(r=u1/(u1+u2)) 
  
  
  # Create grid
  gx <- seq(-1, 1, length = 20)
  gy <- seq(-1, 1, length = 20)


  # Predation
  pdf(paste0(path, "../", "baselineH_pred", ".pdf"), width = 8, height = 8)
  # Interpolate using akima-style linear interpolation
  interp <- akima::interp(df$qm, df$qf, 1-df$pred, xo = gx, yo = gy, duplicate = "mean")
  # Heatmap
  image.plot(
    interp$x, interp$y, interp$z,
    col  = hcl.colors(100, "Inferno"),
    zlim = c(0,1),
    xlab = expression(q[m]),
    ylab = expression(q[f]),
    legend.lab = expression("Nest defence")
  )
  # Data support
  # points(df$qm, df$qf, pch = 16, cex = 0.2, col = rgb(1,1,1,0.3))
  dev.off()

  # Effort
  pdf(paste0(path, "../", "baselineH_u", ".pdf"), width = 8, height = 8)
  # Interpolate using akima-style linear interpolation
  interp <- akima::interp(df$qm, df$qf, df$r, xo = gx, yo = gy, duplicate = "mean")
  # Heatmap
  image.plot(
    interp$x, interp$y, interp$z,
    col  = hcl.colors(100, "Inferno"),
    zlim = c(0,1),
    xlab = expression(q[m]),
    ylab = expression(q[f]),
    legend.lab = expression("male:female effort ratio")
  )
  # Data support
  # points(df$qm, df$qf, pch = 16, cex = 0.2, col = rgb(1,1,1,0.3))
  dev.off()

  file.copy(
    paste0(path, "../", "baselineH_u", ".pdf"),
    paste0(path, "../../", "baselineH_u", ".pdf"),
    overwrite = TRUE
  )
  file.copy(
    paste0(path, "../", "baselineH_pred", ".pdf"),
    paste0(path, "../../", "baselineH_pred", ".pdf"),
    overwrite = TRUE
  )
}

run_trait_plot <- function(
  path,
  x_var,
  x_label,
  out_name,
  y_limits = c(-1, 1)
) {

  readfile <- read.csv(paste0(path, "summaries.csv")) %>%
    mutate(b_s = 1/(1-b_s)) %>%
    filter(b_s <= 100)

  x_sym <- rlang::sym(x_var)
  df_fecundity <- readfile %>%
    mutate(fast_slow = mean_f_s_w, slow_slow = s_w, fast_fast = f_w) %>%
    dplyr::select(!!x_sym, fast_slow, fast_fast, slow_slow) %>%
    tidyr::gather(Loci, Value, -!!x_sym)
  df_effort <- readfile %>%
    mutate(
      fast  = mean_faster_effort,
      slow  = mean_slower_effort
    ) %>%
    select(!!x_sym, fast, slow) %>%
    gather(Loci, Value, -!!x_sym)
  df_info <- readfile %>%
    mutate(
      fast  = mean_fast_h,
      slow  = mean_slow_h
    ) %>%
    select(!!x_sym, fast, slow) %>%
    gather(Loci, Value, -!!x_sym)
  df_sociality <- readfile %>%
    mutate(
      fast  = gamma * (1-mean_fast_h),
      slow  = gamma * (1-mean_slow_h)
    ) %>%
    select(!!x_sym, fast, slow) %>%
    gather(Loci, Value, -!!x_sym)
  
  # Create combined dataframe for legend with all loci
  df_legend <- bind_rows(
    df_fecundity,
    df_effort,
    df_info
  )
  
  # Define all levels and labels
  legend_levels <- c("fast_slow", "fast_fast", "slow_slow", "fast", "slow")
  legend_labels <- c("Fast-Slow", "Fast-Fast", "Slow-Slow", "Fast", "Slow")
  
  # Set factor levels for all dataframes
  df_fecundity$Loci <- factor(df_fecundity$Loci, levels = legend_levels)
  df_effort$Loci <- factor(df_effort$Loci, levels = legend_levels)
  df_info$Loci <- factor(df_info$Loci, levels = legend_levels)
  df_sociality$Loci <- factor(df_sociality$Loci, levels = legend_levels)
  df_legend$Loci <- factor(df_legend$Loci, levels = legend_levels)
  
  x_limits <- range(readfile[[x_var]])
  
  fecundity <- ggplot(df_fecundity, aes(x = !!x_sym, y = Value, color = Loci, group = Loci, fill  = Loci)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
    scale_color_discrete(labels = legend_labels, drop = FALSE) +
    scale_fill_discrete(labels = legend_labels, drop = FALSE) +
    theme_classic(base_size = 12) +
    xlab(x_label) +
    ylab("Relative fecundity") +
    theme(
      panel.grid.major = element_line(colour = "grey80", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey90", linewidth = 0.2),
      legend.position = "none") +
    coord_cartesian(xlim = x_limits, ylim = c(0,2.5))
  
  effort <- ggplot(df_effort, aes(x = !!x_sym, y = Value, color = Loci, group = Loci, fill  = Loci)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
    scale_color_discrete(labels = legend_labels, drop = FALSE) +
    scale_fill_discrete(labels = legend_labels, drop = FALSE) +
    theme_classic(base_size = 12) +
    xlab(x_label) +
    ylab("Provisioning effort") +
    theme(
      panel.grid.major = element_line(colour = "grey80", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey90", linewidth = 0.2),
      legend.position = "none") +
    coord_cartesian(xlim = x_limits, ylim = c(0,1))
  
  info <- ggplot(df_info, aes(x = !!x_sym, y = Value, color = Loci, group = Loci, fill  = Loci)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
    scale_color_discrete(labels = legend_labels, drop = FALSE) +
    scale_fill_discrete(labels = legend_labels, drop = FALSE) +
    theme_classic(base_size = 12) +
    xlab(x_label) +
    ylab("Relative uncertainty") +
    theme(
      panel.grid.major = element_line(colour = "grey80", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey90", linewidth = 0.2),
      legend.position = "none") +
    coord_cartesian(xlim = x_limits, ylim = c(0,1))
  
  sociality <- ggplot(df_sociality, aes(x = !!x_sym, y = Value, color = Loci, group = Loci, fill  = Loci)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
    scale_color_discrete(labels = legend_labels, drop = FALSE) +
    scale_fill_discrete(labels = legend_labels, drop = FALSE) +
    theme_classic(base_size = 12) +
    xlab(x_label) +
    ylab("Exploitation-synergy continuum") +
    theme(
      panel.grid.major = element_line(colour = "grey80", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey90", linewidth = 0.2),
      legend.position = "none") +
    coord_cartesian(xlim = x_limits, ylim = c(-1,1))
  
  # Create dummy plot for legend extraction
  library(cowplot)
  legend_plot <- ggplot(df_legend, aes(x = !!x_sym, y = Value, color = Loci, fill = Loci)) +
    geom_point() +
    scale_color_discrete(labels = legend_labels, drop = FALSE) +
    scale_fill_discrete(labels = legend_labels, drop = FALSE) +
    theme(legend.position = "bottom",
          legend.title = element_blank())
  
  # Extract legend
  legend <- get_legend(legend_plot)
  
  # Combine plots without legend
  plots_combined <- ggarrange(effort, info, fecundity, sociality, 
                              ncol = 2, nrow = 2, legend = "none")
  
  # Add legend at bottom
  final_plot <- ggarrange(plots_combined, legend, 
                          ncol = 1, nrow = 2, 
                          heights = c(10, 1))

  pdf(paste0(path, "../", out_name, ".pdf"), width = 8, height = 6)
  print(final_plot)
  dev.off()

  file.copy(
    paste0(path, "../", out_name, ".pdf"),
    paste0(path, "../../", out_name, ".pdf"),
    overwrite = TRUE
  )
}

run_mehtrait_plot <- function(
  path,
  x_var,
  x_label,
  out_name,
  y_limits = c(-1, 1)
) {

  readfile <- read.csv(paste0(path, "summaries.csv"))

  # tidy evaluation
  x_sym <- rlang::sym(x_var)


  df_long2 <- readfile %>%
    dplyr::select(!!x_sym,
      mean_faster_effort,
      mean_slower_effort,
      mean_fast_h,
      mean_slow_h,
      effort_dif) %>%
    tidyr::gather(Loci, Value, -!!x_sym)

  df_long3 <- readfile %>%
    dplyr::select(!!x_sym, nu, gamma, lambda) %>%
    tidyr::gather(Loci, Value, -!!x_sym)

  x_limits <- range(readfile[[x_var]])

  phen <- ggplot(df_long2, aes(
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
    xlab("") +
    ylab("Behaviours") +
    theme(
      panel.grid.major = element_line(colour = "grey80", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey90", linewidth = 0.2)
    ) +
    scale_color_discrete(name = "Traits", labels = trait_labels) +
    scale_fill_discrete(name = "Traits", labels = trait_labels) +
    coord_cartesian(xlim = x_limits, ylim = c(0,1))

  loc <- ggplot(df_long3, aes(
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
    xlab("") +
    ylab("Loci values") +
    theme(
      panel.grid.major = element_line(colour = "grey80", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey90", linewidth = 0.2)
    ) +
    scale_color_discrete(name = "Key loci", labels = loci_labels) +
    scale_fill_discrete(name = "Key loci", labels = loci_labels) +
    coord_cartesian(xlim = x_limits, ylim = c(-1,1))
  pdf(paste0(path, "../", out_name, ".pdf"), width = 8, height = 6)
  print((loc/phen))
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

run_c_q_plot <- function(path) {
  run_trait_plot(
    path,
    x_var   = "c_q",
    x_label = bquote("Fast POLS mortality cost "(c[q])),
    out_name = "c_q"
  )
}

run_h_plot <- function(path) {
  run_trait_plot(
    path,
    x_var   = "h",
    x_label = bquote("POLS slope of nest defence "(eta)),
    out_name = "h"
  )
}

run_theta_plot <- function(path) {
  run_trait_plot(
    path,
    x_var   = "theta",
    x_label = bquote("POLS slope of mortality "(theta)),
    out_name = "theta"
  )
}

run_sigma_plot <- function(path) {
  run_trait_plot(
    path,
    x_var   = "sigma0",
    x_label = bquote("Standard deviation of POLS "(sigma)),
    out_name = "sigma0"
  )
}

run_sigmacue_plot <- function(path) {
  run_trait_plot(
    path,
    x_var   = "sigma_cue",
    x_label = bquote("Social cue noise"),
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
    x_label = bquote("Baseline lifespan (years)"),
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
