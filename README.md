# Dealing with spatio-temporal autocorrelation

This is a series of scripts that fits increasingly complex models considering spatial and temporal autocorrelation.

# Scripts -- In the `R` direcory

- `01 Initial plot.R`: Produces Figure 1, an initial plot of the study area.

- `02 Read and clean_ data.R`: Reads and cleans the data. Outputs `df.rds`.

- `03_Fit_non-spatial_model.R`: Fits a mixed-effects model without spatial or temporal autocorrelation. Outputs `df_spatial.rds`, `mesh.rds`, `m1` and Figure 2.

- `04_Fit_spatial_model.R`: Adds spatial autocorrelation. Outputs `m2` and Figure 3.

- `05_Fit_spatial_model_with_ar1_temporal_autocorrelation.R`: Adds temporal autocorrelation (AR1). Outputs `m3` and Figure 4.

- `06_Fit_spatial_model_with_random_walk_temporal_autocorrelation.R`: Adds temporal autocorrelation (random walk). Outputs `m4` and Figure 5.

- `07_Compare_model_performance.R`: Compares model performance and outputs final diagnostics: Figure S1 and Figure S2

- `Helpers.R`: Generic helper functions used in other scripts
