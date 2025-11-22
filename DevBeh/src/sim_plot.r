#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(patchwork)
  library(viridis)
})
readfile <- read.csv(path)

suppressMessages({
  suppressWarnings({
    # Benefit of strategic diversity
    pdf(paste0(path, ".pdf"), width = 8, height = 6)

      # Prepare first dataset
    df_long1 <- readfile %>%
      select(i, u_base, rho, nu, 
            gamma, lambda, c, m,x) %>%
      gather(`Loci`, Value, -i)

    # # Prepare second dataset
    # df_long2 <- readfile %>%
    #   select(i, m, x) %>%
    #   gather(`Loci`, Value, -i)

      # Determine shared x-axis limits
      x_limits <- range(readfile$i)
    # First plot
    p1 <- ggplot(df_long1,
                  aes(x = i, y = Value,
                      color = `Loci`,
                      group = `Loci`)) +
      geom_point(alpha = 0.3, size = 1) +
      geom_smooth(se = FALSE, method = "gam") + # smooths each
      theme_classic(base_size = 12) +
      xlab("Generations") +
      coord_cartesian(xlim = x_limits, ylim = c(-1, 1)) +
      scale_color_viridis_d(option = "turbo")  # different color set

    # Second plot
    # p2 <- ggplot(df_long2,
    #               aes(x = i, y = Value,
    #                   color = `Loci`,
    #                   group = `Loci`)) +
    #   geom_point(alpha = 0.3, size = 1) +
    #   geom_smooth(se = FALSE, method = "gam") + # smooths each
    #   theme_classic(base_size = 12) +
    #   xlab("Generations") +
    #   coord_cartesian(xlim = x_limits) +
    #   scale_color_brewer(palette = "Set2") + 
    #   theme(legend.title = element_blank())

    # Combine plots
    # invisible(print(p1 / p2 + plot_layout(ncol = 1, heights = c(1,1))))
    invisible(print(p1))
    invisible(dev.off())
  })
})