suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
library(patchwork)
})
run_b_s_plot <- function(path) {

  # Read CSV
  readfile <- read.csv(paste0(path, "summaries.csv"))

  # Prepare datasets
  df_long1 <- readfile %>%
    select(b_s, u_base, rho, nu,
           gamma, lambda, c, m,x) %>%
    gather(`Negotiation Loci`, Value, -b_s)

  # df_long2 <- readfile %>%
  #   select(b_s, m,x) %>%
  #   gather(`Observation Loci`, Value, -b_s)

  # Determine shared x-axis limits
  x_limits <- range(readfile$b_s)

  # First plot
  p1 <- ggplot(df_long1, aes(x = b_s, y = Value,
                             color = `Negotiation Loci`,
                             group = `Negotiation Loci`)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "gam") +
    theme_classic(base_size = 12) +
    xlab("Survival benefits of parental care (b_s)") +
    coord_cartesian(xlim = x_limits, ylim = c(-1, 1)) +
    scale_color_viridis_d(option = "turbo")  # different color set

  # Second plot
  # p2 <- ggplot(df_long2, aes(x = b_s, y = Value,
  #                            color = `Observation Loci`,
  #                            group = `Observation Loci`)) +
  #   geom_point(alpha = 0.3, size = 1) +
  #   geom_smooth(se=FALSE, method = "gam", formula = y ~ s(x, bs = "cs")) +
  #   theme_classic(base_size = 12) +
  #   xlab("Survival benefits of parental care (b_s)") +
  #   coord_cartesian(xlim = x_limits) +
  #   scale_color_brewer(palette = "Set2") + # different color set
  #   theme(legend.title = element_blank())

  # Open PDF
  pdf(paste0(path, "b_s.pdf"), width = 8, height = 6)

  # Stack plots with patchwork (same width)
  # print(p1 / p2 + plot_layout(ncol = 1, heights = c(1,1)))
  print(p1)

  # Close PDF
  dev.off()
}

run_b_f_plot <- function(path) {

  # Read CSV
  readfile <- read.csv(paste0(path, "summaries.csv"))

  # Prepare datasets
  df_long1 <- readfile %>%
    select(b_f, u_base, rho, nu,
           gamma, lambda, c, m,x) %>%
    gather(`Loci`, Value, -b_f)

  # Determine shared x-axis limits
  x_limits <- range(readfile$b_f)

  # First plot
  p1 <- ggplot(df_long1, aes(x = b_f, y = Value,
                             color = `Loci`,
                             group = `Loci`)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "gam") +
    theme_classic(base_size = 12) +
    xlab("Benefits of parental care (b_f)") +
    coord_cartesian(xlim = x_limits, ylim = c(-1, 1)) +
    scale_color_viridis_d(option = "turbo")  # different color set

  # Open PDF
  pdf(paste0(path, "b_f.pdf"), width = 8, height = 6)

  # Stack plots with patchwork (same width)
  # print(p1 / p2 + plot_layout(ncol = 1, heights = c(1,1)))
  print(p1)

  # Close PDF
  dev.off()
}

run_sigma_plot <- function(path) {

  # Read CSV
  readfile <- read.csv(paste0(path, "summaries.csv"))

  # Prepare datasets
  df_long1 <- readfile %>%
    select(sigma0, u_base, rho, nu,
           gamma, lambda, c, m,x) %>%
    gather(`Negotiation Loci`, Value, -sigma0)

  # df_long2 <- readfile %>%
  #   select(sigma0, m,x) %>%
  #   gather(`Observation Loci`, Value, -sigma0)

  # Determine shared x-axis limits
  x_limits <- range(readfile$sigma0)

  # First plot
  p1 <- ggplot(df_long1, aes(x = sigma0, y = Value,
                             color = `Negotiation Loci`,
                             group = `Negotiation Loci`)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "gam") +
    theme_classic(base_size = 12) +
    xlab("Standard deviation of pace-of-life phenotype (sigma)") +
    coord_cartesian(xlim = x_limits, ylim = c(-1, 1)) +
    scale_color_viridis_d(option = "turbo")  # different color set

  # Second plot
  # p2 <- ggplot(df_long2, aes(x = sigma0, y = Value,
  #                            color = `Observation Loci`,
  #                            group = `Observation Loci`)) +
  #   geom_point(alpha = 0.3, size = 1) +
  #   geom_smooth(se=FALSE, method = "gam", formula = y ~ s(x, bs = "cs")) +
  #   theme_classic(base_size = 12) +
  #   xlab("Standard deviation of pace-of-life phenotype (sigma)") +
  #   coord_cartesian(xlim = x_limits) +
  #   scale_color_brewer(palette = "Set2") + 
  #   theme(legend.title = element_blank())

  # Open PDF
  pdf(paste0(path, "sigma0.pdf"), width = 8, height = 6)

  # Stack plots with patchwork (same width)
  # print(p1 / p2 + plot_layout(ncol = 1, heights = c(1,1)))
  print(p1)

  # Close PDF
  dev.off()
}

run_sigmacue_plot <- function(path) {

  # Read CSV
  readfile <- read.csv(paste0(path, "summaries.csv"))

  # Prepare datasets
  df_long1 <- readfile %>%
    select(sigma_cue, u_base, rho, nu,
           gamma, lambda, c, m,x) %>%
    gather(`Negotiation Loci`, Value, -sigma_cue)

  # df_long2 <- readfile %>%
  #   select(sigma_cue, m,x) %>%
  #   gather(`Observation Loci`, Value, -sigma_cue)

  # Determine shared x-axis limits
  x_limits <- range(readfile$sigma_cue)

  # First plot
  p1 <- ggplot(df_long1, aes(x = sigma_cue, y = Value,
                             color = `Negotiation Loci`,
                             group = `Negotiation Loci`)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "gam") +
    theme_classic(base_size = 12) +
    xlab("Standard deviation of social cue") +
    coord_cartesian(xlim = x_limits, ylim = c(-1, 1)) +
    scale_color_viridis_d(option = "turbo")  # different color set

  # Second plot
  # p2 <- ggplot(df_long2, aes(x = sigma_cue, y = Value,
  #                            color = `Observation Loci`,
  #                            group = `Observation Loci`)) +
  #   geom_point(alpha = 0.3, size = 1) +
  #   geom_smooth(se=FALSE, method = "gam", formula = y ~ s(x, bs = "cs")) +
  #   theme_classic(base_size = 12) +
  #   xlab("Standard deviation of social cue") +
  #   coord_cartesian(xlim = x_limits) +
  #   scale_color_brewer(palette = "Set2") + 
  #   theme(legend.title = element_blank())

  # Open PDF
  pdf(paste0(path, "sigma_cue.pdf"), width = 8, height = 6)

  # Stack plots with patchwork (same width)
  # print(p1 / p2 + plot_layout(ncol = 1, heights = c(1,1)))
  print(p1)

  # Close PDF
  dev.off()
}

run_cuecost_plot <- function(path) {

  # Read CSV
  readfile <- read.csv(paste0(path, "summaries.csv"))

  # Prepare datasets
  df_long1 <- readfile %>%
    select(c_v, u_base, rho, nu,
           gamma, lambda, c, m,x) %>%
    gather(`Loci`, Value, -c_v)
    
    
  # First plot
  p1 <- ggplot(df_long1, aes(x = c_v, y = Value,
                             color = `Loci`,
                             group = `Loci`)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "gam") +
    theme_classic(base_size = 12) +
    xlab("Cost of social cue (c_v)") +
    coord_cartesian(xlim = c(0, 1), ylim = c(-1, 1)) +
    scale_color_viridis_d(option = "turbo")  # different color set

  # Open PDF
  pdf(paste0(path, "cue_cost.pdf"), width = 8, height = 6)

  # Stack plots with patchwork (same width)
  # print(p1 / p2 + plot_layout(ncol = 1, heights = c(1,1)))
  print(p1)

  # Close PDF
  dev.off()
}

run_mort_plot <- function(path) {

  # Read CSV
  readfile <- read.csv(paste0(path, "summaries.csv"))

  # Prepare datasets
  df_long1 <- readfile %>%
    select(s_p, u_base, rho, nu,
           gamma, lambda, c, m,x) %>%
    gather(`Loci`, Value, -s_p)

  # df_long2 <- readfile %>%
  #   select(s_p) %>%
  #   gather(`Loci`, Value, -s_p)

  # Determine shared x-axis limits
  x_limits <- range(readfile$s_p)

  # First plot
  p1 <- ggplot(df_long1, aes(x = s_p, y = Value,
                             color = `Loci`,
                             group = `Loci`)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "gam") +
    theme_classic(base_size = 12) +
    xlab("Baseline survival rate (s_p)") +
    coord_cartesian(xlim = x_limits, ylim = c(-1, 1)) +
    scale_color_viridis_d(option = "turbo")  # different color set

  # # Second plot
  # p2 <- ggplot(df_long2, aes(x = s_p, y = Value,
  #                            color = `Loci`,
  #                            group = `Loci`)) +
  #   geom_point(alpha = 0.3, size = 1) +
  #   geom_smooth(se=FALSE, method = "gam", formula = y ~ s(x, bs = "cs")) +
  #   theme_classic(base_size = 12) +
  #   xlab("Baseline survival rate (s_p)") +
  #   coord_cartesian(xlim = x_limits) +
  #   scale_color_brewer(palette = "Set2") + 
  #   theme(legend.title = element_blank())

  # Open PDF
  pdf(paste0(path, "s_p.pdf"), width = 8, height = 6)

  # Stack plots with patchwork (same width)
  # print(p1 / p2 + plot_layout(ncol = 1, heights = c(1,1)))
  print(p1)

  # Close PDF
  dev.off()
}

run_divorce_plot <- function(path) {

  # Read CSV
  readfile <- read.csv(paste0(path, "summaries.csv"))

  # Prepare datasets
  df_long1 <- readfile %>%
    select(div_rate, u_base, rho, nu,
           gamma, lambda, c, m,x) %>%
    gather(`Loci`, Value, -div_rate)

  # Determine shared x-axis limits
  x_limits <- range(readfile$div_rate)

  # First plot
  p1 <- ggplot(df_long1, aes(x = div_rate, y = Value,
                             color = `Loci`,
                             group = `Loci`)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "gam") +
    theme_classic(base_size = 12) +
    xlab("Divorce rate") +
    coord_cartesian(xlim = x_limits, ylim = c(-1, 1)) +
    scale_color_viridis_d(option = "turbo")  # different color set


  # Open PDF
  pdf(paste0(path, "div_rate.pdf"), width = 8, height = 6)

  # Stack plots with patchwork (same width)
  # print(p1 / p2 + plot_layout(ncol = 1, heights = c(1,1)))
  print(p1)

  # Close PDF
  dev.off()
}
