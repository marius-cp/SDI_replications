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
  
  These data sets are stored in the subfolder `data`.  
  Currently, it is unclear whether these data sets can be made publicly available due to data licensing restrictions.  
  This script also includes the code used to produce the results for the motivational example.

- **`evaluation_vola.R`**  
  Produces the results for the volatility application.

- **`evaluation_var.R`**  
  Produces the results for the Value-at-Risk application.

- **`timeseriesplot.R`**  
  Generates the time series plot shown in the Appendix.

- **`MZdiagnostic.R`**  
  Produces the Mincer--Zarnowitz regression diagnostic plots shown in the Appendix.

### Notes on data availability
Unfortunately, the data cannot be shared due to licensing restrictions.
