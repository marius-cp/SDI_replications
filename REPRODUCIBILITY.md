# Reproducibility guide

This material reproduces the results in *Statistical Inference for Score Decompositions*. The reference software environment is R 4.2.1 on macOS (Apple silicon); exact package versions are recorded in `SESSION_INFO.txt`. `00_functions/check_dependencies.R` reports version differences but allows the replication to continue. It stops only when a required package is missing.

The reference run used `SDI` version 0.0.1 from commit `0edf793ea183061d3d6eb53cbb60faaffb5b23e9` of `marius-cp/SDI`. Other installed versions are reported but do not stop the replication attempt. Scripts never install or update packages automatically.

On a clean R 4.2.1 installation, run `Rscript install_dependencies.R` once to install the recorded package versions. This setup script is never invoked by the analysis scripts.

## Recommended run order

`master.R` lists every replication script in execution order. Select the desired parts with the `TRUE`/`FALSE` switches at the top of that file, then use **Source** in RStudio or run:

```sh
Rscript master.R
```

The supplied settings rebuild the replication from the available raw inputs and regenerate the paper figures and printed results. The three simulation scripts took approximately 15 hours in total in the reference run.

Keep `RUN_MAIN_SIMULATIONS <- TRUE` to rerun all Monte Carlo simulations. In the reference run, the mean simulation took approximately 6 hours, the main quantile simulation approximately 4.7 hours, and the additional, unreported finite-sample exercise approximately 4.5 hours. The main quantile simulation writes four stable files named `sim_q_parameterized_part1.rds` through `sim_q_parameterized_part4.rds`; no date is added to filenames and the temporary combined object is not written to disk.

Keep `RUN_APPLICATION_PREPARATION <- TRUE` for a rebuild of the application inputs. If the licensed Refinitiv input described in `DATA_AVAILABILITY.md` is absent, `master.R` reports this and skips only `cleaning.R` and `forecasting.R`. Financial figures and tables are still regenerated when the prepared financial RDS files are available; otherwise that output section is also skipped with a message.

The other switches separately control the simulation figures and the application figures and printed results. Each individual analysis script can also still be run on its own.

The additional, unreported finite-sample quantile exercise is included when `RUN_MAIN_SIMULATIONS <- TRUE`. It can also be run separately with `Rscript 01_simulation/sim_quantile_unreported.R`.

## Paper-output crosswalk

| Paper output | Replication script | Main input |
|---|---|---|
| `sim_m.pdf` | `01_simulation/sim_mean_plots.R` | `sim_m_parameterized.rds` |
| `sim_q.pdf` | `01_simulation/sim_quantile_plots.R` | four `sim_q_parameterized_part*.rds` files |
| `m_true.pdf`, `q_true.pdf`, `q_xi0_choice.pdf` | `01_simulation/plot_true_components.R` | closed-form expressions and `sim_q_compare.rds` |
| `appl_inflation_timeseries.pdf`, `appl_infl.pdf` | `02_application/inflation/inflation.R` | `inflation_mean.rds` |
| `MZdiagnostics_inflation.pdf` | `02_application/inflation/MZdiagnostic.R` | `inflation_mean.rds` |
| `appl_vola_emini_comb.pdf` | `02_application/finrisk/evaluation_vola.R` | `emini_fcasts_vola.rds` |
| `table_1_upper_variance.tex` (upper part of manuscript Table 1) | `02_application/finrisk/evaluation_vola.R` | `emini_fcasts_vola.rds` |
| `appl_VaR_emini_comb.pdf` | `02_application/finrisk/evaluation_var.R` | two `emini_fcasts_VaR_*.rds` files |
| `table_1_lower_var.tex` (lower part of manuscript Table 1) | `02_application/finrisk/evaluation_var.R` | `emini_fcasts_VaR_1.rds` |
| `table_5_var_backtests.tex` (manuscript Table 5) | `02_application/finrisk/evaluation_var.R` | two `emini_fcasts_VaR_*.rds` files |
| `timeseries.pdf` | `02_application/finrisk/timeseriesplot.R` | prepared financial forecast files |
| `MZdiagnostics_vola.pdf`, `MZdiagnostics_VaR.pdf` | `02_application/finrisk/MZdiagnostic.R` | prepared financial forecast files |

The inflation p-values printed by the reference run are 0.811 (equal MCB), 0.037 (equal DSC), and 0.284 (equal overall score), matching the manuscript. The remaining application tables are printed by the two evaluation scripts in the same run.
