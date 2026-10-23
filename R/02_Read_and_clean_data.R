# Read the data for the analysis, and sort it out
  # Do any data manipulation here
    # Written by Dave S
      # October 2024


# Attach package libraries ------------------------------------------------

  pacman::p_load(tidyverse)


# Get the data ------------------------------------------------------------

 df <- read_csv("Data/blob_kelp.csv") %>% # Read in the data
    mutate(across(where(is.character), factor), # Make character variables into factors
          pixelid = factor(pixelid), # Convert this one manually, because it is a number, not a character
          region = fct_relevel(region, "north"), # Make "north" the first factor level
          period = fct_relevel(period, "before")) # Make "before" the first factor level

 
# Save the data -----------------------------------------------------------

   write_rds(df, "output/df.rds") # Every script ends by producing outputs that become inputs for the next script

 