# Compare model performance and check final diagnostics
  # Written by Dave S
    # October 2024


# Attach package libraries ------------------------------------------------

pacman::p_load(
  tidyverse,
  performance
  )
  source("R/Helpers.R")


# Compare models ----------------------------------------------------------

  m1 <- read_rds("output/m1.rds")
  m2 <- read_rds("output/m2.rds")
  m3 <- read_rds("output/m3.rds")
  m4 <- read_rds("output/m4.rds")
  
  compare_performance(m1, m2, m3, m4)
  

# Partial autocorrelation plots for m4 ------------------------------------

  plot_pacf_by_group(model = m4, 
                     data = dat, 
                     time_var = year, 
                     grouping_var = pixelid)
  ggsave("figures/Fig_S1.pdf")


# DHARMA residuals --------------------------------------------------------
  
  sdmTMB_residuals(m4) # This uses a helper function
  ggsave("figures/Fig_S2.pdf")
  