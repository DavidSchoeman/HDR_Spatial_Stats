# Fit a modelaccounting for spatial and temporal autocorrelation (AR1)
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


# Build the spatio-temporal model (AR1) -----------------------------------

  m3 <- sdmTMB(area ~ region * period +
              (1 | pixelid),
             mesh = mesh,
             spatial = "on",
             time = "year", # Specify that year controls for time
             spatiotemporal = "ar1", # Specify autoregressive correlation of order 1
             data = dat,
             family = tweedie())
  tidy(m3) # Tidy model summary
  Anova(m3) # Rough check for significance

  
# Check diagnostics -------------------------------------------------------

  sdmTMB_residuals(m3) # This uses a helper function
  

# Plot outputs ------------------------------------------------------------

  # Predict values from our model for this prediction tibble
    preds3 <- estimate_expectation(m3, new_dat) %>% 
      tibble()

  # Plot the results
    ggplot() +
      geom_errorbar(data = preds3,
                    aes(x = period, ymin = CI_low, ymax = CI_high, colour = region),
                    width = .05,
                    linewidth = 1,
                    position = position_dodge(.3)) +
      geom_point(data = preds3,
                 aes(x = period, y = Predicted, colour = region),
                 position = position_dodge(.3),
                 size = 3) +
      labs(colour = "Region",
           y = "Areal coverage of kelp ± 95% confidence (ha)",
           x = "Timing relative to the onset of MHWs") +
      scale_colour_manual(values = col_pal)

    ggsave("figures/Fig_4.pdf") # Output plot
    write_rds(m3, "output/m3.rds") # Keep a copy of the model
