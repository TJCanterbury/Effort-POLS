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

    invisible(png(paste0(path, ".png"), width = 720, height = 480))

    # Prepare first dataset
    df_long1 <- readfile %>%
      select(i, u_base, rho, nu, 
            gamma, lambda, c, m,x) %>%
      gather(`Loci`, Value, -i)

      x_limits <- range(readfile$i)
      
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

    invisible(print(p1))
    invisible(dev.off())
  })
})