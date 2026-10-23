# Fit a model without consideration of spatial or temporal autocorrelation
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

  dat <- read_rds("output/df.rds") %>% 
    add_utm_columns(ll_names = c("lon", "lat")) # Add columns of UTM coordinates...we'll need this soon
  write_rds(dat, "output/df_spatial.rds") # Save the spatial tibble


# Make SPDE mesh ----------------------------------------------------------

  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 10)
  plot(mesh)
  write_rds(mesh, "output/mesh.rds") # Store this for later use
  

# Build the non-spatial model ---------------------------------------------

  m1 <- sdmTMB(area ~ region * period +
                 (1 | pixelid),
               mesh = mesh,
               spatial = "off",
               data = dat,
               family = tweedie())
  
  tidy(m1) # Tidy model summary
  Anova(m1) # Rough check for significance

  
# Check diagnostics -------------------------------------------------------

  sdmTMB_residuals(m1) # This uses a helper function
  

# Plot outputs ------------------------------------------------------------
  
  # Build a prediction tibble
    new_dat <- expand_grid(region = levels(dat$region),
                           period = levels(dat$period), # A tibble with two variables containing all possible combinations of region and period 
                           year = min(dat$year)) # We'll add a single value for year, here, even though we're not using it yet
  
    write_rds(new_dat, "output/new_dat.rds") # Take a copy

  # Predict values from our model for this prediction tibble
    preds1 <- estimate_expectation(m1, new_dat) %>% 
      tibble()

  # Plot the results
    ggplot() +
      geom_errorbar(data = preds1,
                    aes(x = period, ymin = CI_low, ymax = CI_high, colour = region),
                    width = .05,
                    linewidth = 1,
                    position = position_dodge(.3)) +
      geom_point(data = preds1,
                 aes(x = period, y = Predicted, colour = region),
                 position = position_dodge(.3),
                 size = 3) +
      labs(colour = "Region",
           y = "Areal coverage of kelp ± 95% confidence (ha)",
           x = "Timing relative to the onset of MHWs") +
      scale_colour_manual(values = col_pal)

    ggsave("figures/Fig_2.pdf") # Output plot
    write_rds(m1, "output/m1.rds") # Keep a copy of the model
