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

  dat <- read_rds("data/df.rds") %>% 
    add_utm_columns(ll_names = c("lon", "lat")) # Add columns of UTM coordinates...we'll need this soon
  write_rds(dat, "data/df_spatial.rds") # Save the spatial tibble


# Make SPDE mesh ----------------------------------------------------------

  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 10)
  plot(mesh)
  write_rds(mesh, "data/mesh.rds") # Store this for later use
  

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
  ggsave("figures/Fig_2.pdf")


# Deploy the function


m2 <- sdmTMB(area ~ region * period,
             mesh = mesh,
             spatial = "on", # Switch on the random spatial mesh
             data = dat,
             family = tweedie())
Anova(m2)

# Predict values from our model for this prediction tibble
  preds2 <- estimate_expectation(m2, new_dat) %>% 
    tibble()

# Plot the results
  ggplot(preds2) +
    geom_errorbar(aes(x = period, ymin = CI_low, ymax = CI_high, colour = region),
                  width = .05,
                  linewidth = 1,
                  position = position_dodge(.3)) +
    geom_point(aes(x = period, y = Predicted, colour = region),
               position = position_dodge(.3),
               size = 3) +
    labs(colour = "Region",
         y = "Areal coverage of kelp ± 95% confidence (ha)",
         x = "Timing relative to the onset of MHWs") +
    scale_colour_manual(values = col_pal)

sdmTMB_residuals(m2)

# Run a Durbin-Watson test for each time series, each time extracting just the p-value from the test
dat %>% 
  arrange(year) %>%  # Sort by year to ensure that the data are ordered in time
  group_by(pixelid) %>% # Group by pixel to ensure that we work one time series at a time
  summarise(p = dwtest(area ~ 1) %>% pluck("p.value")) %>% # For each time series, run the test and extract the p-value
  mutate(significance = ifelse(p <= 0.05, "***significant***", "non-significant")) %>% # Make a new column to make it easier to spot significant results
  print(n = nrow(.)) # Print all rows

m3 <- sdmTMB(area ~ region * period +
             (1 | pixelid),
             mesh = mesh,
             spatial = "on",
             time = "year", # Specify that year controls for time
             spatiotemporal = "ar1", # Specify autoregressive correlation of order 1
             data = dat,
             family = tweedie())
Anova(m3)

sdmTMB_residuals(m3)

# Predict values from our model for this prediction tibble
  preds3 <- estimate_expectation(m3, new_dat) %>% 
    tibble()

# Plot the results
  ggplot(preds3) +
    geom_errorbar(aes(x = period, ymin = CI_low, ymax = CI_high, colour = region),
                  width = .05,
                  linewidth = 1,
                  position = position_dodge(.3)) +
    geom_point(aes(x = period, y = Predicted, colour = region),
               position = position_dodge(.3),
               size = 3) +
    labs(colour = "Region",
         y = "Areal coverage of kelp ± 95% confidence (ha)",
         x = "Timing relative to the onset of MHWs") +
    scale_colour_manual(values = col_pal)

m4 <- sdmTMB(area ~ region * period +
               (1 | pixelid),
             mesh = mesh,
             spatial = "on",
             time = "year", # Specify that year controls for time
             spatiotemporal = "rw", # Specify random walk for temporal autocorrelation
             data = dat,
             family = tweedie())
Anova(m4)

# Predict values from our model for this prediction tibble
  preds4 <- estimate_expectation(m4, new_dat) %>% 
    tibble()

# Plot the results
  ggplot(preds4) +
    geom_errorbar(aes(x = period, ymin = CI_low, ymax = CI_high, colour = region),
                  width = .05,
                  linewidth = 1,
                  position = position_dodge(.3)) +
    geom_point(aes(x = period, y = Predicted, colour = region),
               position = position_dodge(.3),
               size = 3) +
    labs(colour = "Region",
         y = "Areal coverage of kelp ± 95% confidence (ha)",
         x = "Timing relative to the onset of MHWs") +
    scale_colour_manual(values = col_pal)

sdmTMB_residuals(m4)

compare_performance(m1, m2, m3, m4)

plot_pacf_by_group <- function(model, data, time_var, grouping_var) {
df <- data %>% 
  mutate(residuals = sdmTMB_residuals(model, plot = FALSE) %>% 
           pluck("scaledResiduals")) %>%  # Append the residuals from the model to the data
  arrange({{time_var}}) %>% 
  group_by({{grouping_var}}) %>% # Group by pixelid (i.e., by time series)
  summarise(PACF = list(pacf(residuals, plot = FALSE) %>% # For each time series, extract the pacf stats
                          pluck("acf")), # Retaining ONLY the autocorrelation function values
            n = n()) %>% # Compute the number of observations by group
  unnest(cols = PACF) %>%  # Spread the nested tibbles created from the summarise function
  group_by({{grouping_var}}) %>% 
  mutate(lag = row_number()) %>% # Add lags by group
  ungroup()
ggplot() +
  geom_col(data = df, # Column graphs
           aes(x = lag, y = PACF),  # pacf by lag
           fill = "skyblue") +
  facet_wrap(vars({{grouping_var}})) + # Creates separate panels for each ID
  geom_hline(yintercept = 0) + # Horizontal line at zero (i.e., no autocorrelation)
  geom_hline(yintercept = c(1.96, -1.96)/df %>% pull(n) %>% min() %>% sqrt(), color = "red", linetype = "dashed") + # Approx confidence limits
  labs(x = "Lag", y = "Partial Autocorrelation") +
  scale_x_continuous(breaks = breaks_width(2)) +
  theme(axis.text = element_text(size = 5)) # Make x-axis text small enough to prevent overlap
  }  
  
plot_pacf_by_group(model = m4, 
                   data = dat, 
                   time_var = year, 
                   grouping_var = pixelid)

dat %>% 
  mutate(residuals = sdmTMB_residuals(m4, plot = FALSE) %>% 
           pluck("scaledResiduals")) %>% 
  group_by(pixelid) %>% 
  summarise(p = dwtest(residuals ~ 1) %>% pluck("p.value")) %>%
  mutate(significance_p = ifelse(p <= 0.05, "***significant***", "non-significant")) %>% 
  print(n = nrow(.))

sessionInfo()

knitr::purl("SpatialAuto_Session4.qmd", documentation = 0, output = "CodeFromNotes/SpatialAuto_Session4.R")

unlink(list.files("Output", full.names = TRUE, recursive = TRUE), recursive = TRUE)
