# Fit a modelaccounting for spatial and temporal autocorrelation (random walk)
    # Written by Dave S
      # October 2024


# Attach package libraries ------------------------------------------------

  pacman::p_load(
    tidyverse,
    DHARMa,
    sdmTMB,
    lmtest,
    performance,
    modelbased,
    car,
    scales
  )
  source("R/Helpers.R")


# Get the data ------------------------------------------------------------

  dat <- read_rds("output/df_spatial.rds")
  mesh <- read_rds("output/mesh.rds")
  new_dat <- read_rds("output/new_dat.rds")


# Build the spatio-temporal model (rw) -----------------------------------

  m4 <- sdmTMB(area ~ region * period +
             (1 | pixelid),
             mesh = mesh,
             spatial = "on",
             time = "year", # Specify that year controls for time
             spatiotemporal = "rw", # Specify random walk for temporal autocorrelation
             data = dat,
             family = tweedie())
  tidy(m4) # Tidy model summary
  Anova(m4) # Rough check for significance

  
# Check diagnostics -------------------------------------------------------

  sdmTMB_residuals(m4) # This uses a helper function
  

# Plot outputs ------------------------------------------------------------

  # Predict values from our model for this prediction tibble
    preds4 <- estimate_expectation(m4, new_dat) %>% 
      tibble()

  # Plot the results
    ggplot() +
      geom_errorbar(data = preds4,
                    aes(x = period, ymin = CI_low, ymax = CI_high, colour = region),
                    width = .05,
                    linewidth = 1,
                    position = position_dodge(.3)) +
      geom_point(data = preds4,
                 aes(x = period, y = Predicted, colour = region),
                 position = position_dodge(.3),
                 size = 3) +
      labs(colour = "Region",
           y = "Areal coverage of kelp ± 95% confidence (ha)",
           x = "Timing relative to the onset of MHWs") +
      scale_colour_manual(values = col_pal)

    ggsave("figures/Fig_5.pdf") # Output plot
    write_rds(m4, "output/m4.rds") # Keep a copy of the model
