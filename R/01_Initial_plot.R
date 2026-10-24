# Read the data and create an initial plot
    # Written by Dave S
      # October 2024


# Attach package libraries ------------------------------------------------

  pacman::p_load(
    tidyverse,
    terra,
    sf, 
    rnaturalearth
    )


# Get the data to plot ---------------------------------------------------

  df <- read_csv("Data/blob_kelp.csv") %>% # Read in the data
    filter(year == 2006) %>% # Get just the first year, so that we don't have replicated 0.25º grid squares
    dplyr::select(pixelid, lon, lat, region)  # Get just the coordinates of the 0.25º grid squares, and their identities


# Produce the plot ----------------------------------------------------------

  # Make a polygon
   sf::sf_use_s2(FALSE) # Switches off the use of the Google S2 spherical geometry library within the sf package and reverts to using the GEOS library for geometric operations. 
  bbox <- ext(-118, -112, 26.75, 33) # A rough bounding box for the study area
  land <- ne_countries(scale = 50) %>% # A medium-resolution database of country polygons
    st_crop(bbox) %>% # Cropped to our bounding box
    st_geometry() # Retain only the geometries, not any of the other descriptors

  # Plot    
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

  # Save
    ggsave("figures/Fig_1.pdf") # Output to logical directory
    