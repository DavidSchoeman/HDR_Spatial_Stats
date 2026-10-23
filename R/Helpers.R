# Helper functions
  # Written by Dave S
    # October 2024



# Colour palette for plots ------------------------------------------------

  col_pal <- c(c("cadetblue4", "darkorange")) 


# A function for plotting residuals from sdmTMB models --------------------

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
  
  
  
# A function for partial autocorrelation plots ----------------------------

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