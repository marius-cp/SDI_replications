# Statistical Inference for Score Decompositions

Replication material for:

Timo Dimitriadis and Marius Puke (2026), *Statistical Inference for Score Decompositions*. Preprint available on [arXiv](https://arxiv.org/abs/2603.04275).

## Contents

- `00_functions`: shared simulation, plotting, and backtesting functions.
- `01_simulation`: mean and quantile Monte Carlo experiments and the scripts that generate the paper figures.
- `02_application/inflation`: U.S. CPI inflation application and source data preparation.
- `02_application/finrisk`: volatility and Value-at-Risk application. The required Refinitiv data are license-restricted; see `DATA_AVAILABILITY.md`.
- `master.R`: master file that runs the individual replication scripts in their required order.

## Quick start

The reference run used R 4.2.1 and the package versions in `SESSION_INFO.txt`. In particular, it used `SDI` 0.0.1 at GitHub commit `0edf793ea183061d3d6eb53cbb60faaffb5b23e9`. Different installed versions are reported but do not stop the replication attempt. Missing packages still cause an early, informative stop. The analysis scripts do not install or update packages.

From a clean R 4.2.1 installation, install the reference package versions explicitly with `Rscript install_dependencies.R`.

Open `master.R` and select the desired parts using the `TRUE`/`FALSE` switches at the top. The supplied settings rebuild the replication from the available raw inputs and regenerate the paper outputs. If the licensed Refinitiv input is absent, the master file reports this and skips only the financial-data preparation. Financial figures and tables are still regenerated when the prepared financial RDS files are available. The three simulation scripts took approximately 15 hours in total in the reference run. Run the master file in RStudio with **Source**, or from a terminal with:

```sh
Rscript master.R
```

See `REPRODUCIBILITY.md` for the expensive full-simulation workflow and `DATA_AVAILABILITY.md` for the data-access limitation.

## Output safety

Scripts save figures only inside this replication folder. Dropbox/Overleaf export calls are disabled and copying a figure into the paper project is always a manual action.
