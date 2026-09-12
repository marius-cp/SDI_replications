# Data availability and file manifest

## Simulations

The reference simulation output in `01_simulation/data` is included in the repository and can be regenerated entirely by the supplied R code. The large quantile result is split into four files during the run.

## Inflation application

The source files are included in `02_application/inflation/murphy_replication`. The master run creates `inflation_mean.rds` from these inputs. The preparation script records the original sources in its comments, and `02_application/inflation/README.md` gives the literature reference.

## Financial-risk application

The raw one-minute E-mini S&P 500 futures data were obtained from Refinitiv and cannot be redistributed under the current license. The local `02_application/finrisk/data` directory therefore remains excluded from Git.

For an authorized replication, place the raw file at:

`02_application/finrisk/data/index_futures_ESc1_1min_2000_2022.rds`

Then run `cleaning.R` followed by `forecasting.R`, or set `RUN_APPLICATION_PREPARATION <- TRUE` in `master.R` and run the master file.

The reference local files have the following dimensions:

| File | Rows | Columns |
|---|---:|---:|
| `emini_clean.rds` | 5,738 | 13 |
| `emini_fcasts_vola.rds` | 3,617 | 8 |
| `emini_fcasts_VaR_1.rds` | 3,617 | 9 |
| `emini_fcasts_VaR_5.rds` | 3,617 | 9 |
