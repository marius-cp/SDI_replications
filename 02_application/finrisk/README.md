## Replication Material for Application: Financial Volatility and Value-at-Risk Forecast Performance

This subfolder contains the R scripts used for the empirical application on financial volatility and Value-at-Risk (VaR) forecast performance.

### Overview of scripts

- **`cleaning.R`**  
  Cleans the raw Refinitiv data (which cannot be made publicly available) and constructs the return series as well as realized measures.

- **`forecasting.R`**  
  Implements the forecasting models for the volatility and VaR applications and creates the following data sets:
  
  - `emini_fcasts_vola.rds` for the volatility application  
  - `emini_fcasts_VaR_1.rds` for the 1% VaR application  
  - `emini_fcasts_VaR_5.rds` for the 5% VaR application  
  
  These data sets are stored locally in the subfolder `data` and remain outside Git because they are derived from the license-restricted Refinitiv input.
  This script also includes the code used to produce the results for the motivational example.

- **`evaluation_vola.R`**  
  Produces the results for the volatility application and writes the upper, variance-forecast part of manuscript Table 1 to `tables/table_1_upper_variance.tex`.

- **`evaluation_var.R`**  
  Produces the results for the Value-at-Risk application, writes the lower, VaR-forecast part of manuscript Table 1 to `tables/table_1_lower_var.tex`, and writes manuscript Table 5 to `tables/table_5_var_backtests.tex`.

- **`timeseriesplot.R`**  
  Generates the time series plot shown in the Appendix.

- **`MZdiagnostic.R`**  
  Produces the Mincer--Zarnowitz regression diagnostic plots shown in the Appendix.

### Notes on data availability

The Refinitiv input cannot be shared under the current license. The raw and prepared financial data remain outside Git. Exact expected filenames, dimensions, and checksums are recorded in `../../DATA_AVAILABILITY.md` so an authorized replicator can validate a local copy.

The analysis scripts do not install packages and do not write to Dropbox or Overleaf. In `master.R`, keep `RUN_APPLICATION_OUTPUTS <- TRUE` to regenerate all application outputs from prepared local data.
