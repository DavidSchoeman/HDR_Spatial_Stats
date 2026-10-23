library(tidyverse)
library(terra)
library(sf)
library(DHARMa)
library(sdmTMB)
library(lmtest)
library(performance)
library(modelbased)
library(rnaturalearth)
library(car)
library(scales)

# library(tidyverse)
# library(terra)
# library(sf)
# library(DHARMa)
# library(sdmTMB)
# library(lmtest)
# library(performance)
# library(modelbased)
# library(rnaturalearth)
# library(car)
# library(scales)

# Get the data
 df <- read_csv("Data/blob_kelp.csv") %>% # Read in the data
    filter(year == 2006) %>% # Get just the first year, so that we don't have replicated 0.25º grid squares
    dplyr::select(pixelid, lon, lat, region) # Get just the coordinates of the 0.25º grid squares, and their identities

# Build a land polygon
  sf::sf_use_s2(FALSE) # Switches off the use of the Google S2 spherical geometry library within the sf package and reverts to using the GEOS library for geometric operations. 
  bbox <- ext(-118, -112, 26.75, 33) # A rough bounding box for the study area
  land <- ne_countries(scale = 50) %>% # A medium-resolution database of country polygons
    st_crop(bbox) %>% # Cropped to our bounding box
    st_geometry() # Retain only the geometries, not any of the other descriptors

# Build a plot    
  ggplot() + 
    geom_sf(data = land, fill = "wheat") + # Plot the land and colour it in
    geom_tile(data = df, aes(lon, lat, fill = region), # Add equal-sized tiles for the coordinates
              alpha = .2, # Transparent
              colour = "black", # With black outlines
              show.legend = FALSE) + # No legend
    geom_text(data = df, aes(lon, lat, label = pixelid), # Add pixelid labels
              size = 2) + # Make them readable
    geom_hline(yintercept = 29.72, # A horizontal line separating north and central regions
               colour = "red", # Make it red
               lty = "dashed") + # Make it dashed
    scale_fill_manual(values = c("darkorange", "cadetblue4"), # Colour the grid squares by weatther they are in the northern or southern sub-region
      na.value = "transparent") + # Set NAs to be transparent
    labs(x = "", y = "") + # No x- or y-axis labels
    theme(panel.background = element_rect(fill = "lightblue", color = NA)) + # Make the sea blue
    ylim(c(27.05, 32.7)) + # Clip the figure to avoid polygon bounds being visible
    xlim(c(-117.5, -112.5)) # Clip the figure to avoid polygon bounds being visible

dat <- read_csv("Data/blob_kelp.csv") %>% 
  mutate(across(where(is.character), factor), # Make character variables into factors
         pixelid = factor(pixelid), # Convert this one manually, because it is a number, not a character
         region = fct_relevel(region, "north"), # Make "north" the first factor level
         period = fct_relevel(period, "before")) %>% # Make "before" the first factor level 
  add_utm_columns(ll_names = c("lon", "lat")) # Add columns of UTM coordinates...we'll need this later

col_pal <- c(c("cadetblue4", "darkorange")) # The same colour palette as we used in the map
ggplot(dat) +
  geom_line(aes(x = year, y = area, colour = region)) + # Plot the time series, coloured by region
  geom_vline(xintercept = 2014.5, # A vertical line to indicate the separation of "before" and "after" the onset of the MHWs
             linetype = "dashed",
             colour = "red") +
  labs(colour = "Region", # Sort out axis and legend labels
       x = "",
       y = "Areal coverage of kelp (ha)") +
  scale_colour_manual(values = col_pal) + # Set the colour palette
  facet_wrap(vars(pixelid)) + # Make a separate panel for each time series
  theme(legend.position = "bottom", # Put legend at the bottom
        axis.text.x = element_text(size = 5)) # Make x-axis text small enough to prevent overlap

mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 10)
plot(mesh)

m1 <- sdmTMB(area ~ region * period +
               (1 | pixelid),
             mesh = mesh,
             spatial = "off",
             data = dat,
             family = tweedie())

Anova(m1)

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

# A function for plotting residuals from sdmTMB models
sdmTMB_residuals <- function(model, nsim = 1000, plot = TRUE) {
  set.seed(123456) # Set the random seed so that we get the same results for randomisations
  sims <- simulate({{model}}, nsim = nsim, type = "mle-mvn") # We can't use DHARMa residuals "As is" with sdmTMB
  res <- dharma_residuals({{model}}, simulated_response = sims, integerResponse = TRUE, return_DHARMa = TRUE)
  if(plot != TRUE) {
    return(res)
    } else {
      plot(res)
      }
  }

# Deploy the function
sdmTMB_residuals(m1)

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
