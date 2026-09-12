# Master file for the replication of
# "Statistical Inference for Score Decompositions"

# Select the parts to run --------------------------------------------------

# The settings below rebuild the complete replication from the raw inputs.
RUN_MAIN_SIMULATIONS <- TRUE         # Very time-consuming (about 15 hours)
RUN_APPLICATION_PREPARATION <- TRUE  # Refinitiv preparation is skipped if its raw file is absent
RUN_SIMULATION_OUTPUTS <- TRUE
RUN_APPLICATION_OUTPUTS <- TRUE


# Set up the project -------------------------------------------------------

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  project_root <- dirname(normalizePath(sub("^--file=", "", script_arg)))
} else if (
  requireNamespace("rstudioapi", quietly = TRUE) &&
  rstudioapi::isAvailable()
) {
  project_root <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  project_root <- normalizePath(getwd())
}
setwd(project_root)

source("00_functions/check_dependencies.R")

required_public_inputs <- c(
  "02_application/inflation/murphy_replication/redbk32.xls",
  "02_application/inflation/murphy_replication/Individual_CPI.xlsx",
  "02_application/inflation/murphy_replication/cpiQvMd.xlsx"
)
missing_public_inputs <- required_public_inputs[!file.exists(required_public_inputs)]
if (length(missing_public_inputs)) {
  stop(
    "Missing raw input file(s): ",
    paste(missing_public_inputs, collapse = ", "),
    call. = FALSE
  )
}

financial_raw_input <-
  "02_application/finrisk/data/index_futures_ESc1_1min_2000_2022.rds"
has_financial_raw_input <- file.exists(financial_raw_input)

# Each analysis script starts from a clean workspace. Running it in a separate
# R session preserves that behavior and prevents one script from removing the
# objects and switches in this master file.
run_script <- function(path, env = character()) {
  message("\n==> ", path)
  status <- system2(
    command = file.path(R.home("bin"), "Rscript"),
    args = normalizePath(path),
    env = env
  )
  if (!identical(status, 0L)) {
    stop("Script failed: ", path, call. = FALSE)
  }
}


# 1. Main Monte Carlo simulations -----------------------------------------

if (RUN_MAIN_SIMULATIONS) {
  run_script("01_simulation/sim_mean_parameterized.R") # 6 hours
  run_script("01_simulation/sim_quantile_parameterized.R") # 4.7 hours
  run_script("01_simulation/sim_quantile_unreported.R") # 4.5 hours
  run_script(
    "01_simulation/plot_true_components.R",
    env = "SDI_RECOMPUTE_TRUE_COMPONENTS=true"
  )
}


# 2. Application data preparation ----------------------------------------

if (RUN_APPLICATION_PREPARATION) {
  run_script("02_application/inflation/murphy_replication/prepare_inflation.R")
  if (has_financial_raw_input) {
    run_script("02_application/finrisk/cleaning.R")
    run_script("02_application/finrisk/forecasting.R")
  } else {
    message(
      "\nLicensed Refinitiv input not found; financial data preparation was skipped: ",
      financial_raw_input
    )
  }
}


# 3. Simulation figures ---------------------------------------------------

if (RUN_SIMULATION_OUTPUTS) {
  if (!RUN_MAIN_SIMULATIONS) {
    run_script("01_simulation/plot_true_components.R")
  }
  run_script("01_simulation/sim_mean_plots.R")
  run_script("01_simulation/sim_quantile_plots.R")
}


# 4. Application figures and printed results ------------------------------

if (RUN_APPLICATION_OUTPUTS) {
  run_script("02_application/inflation/inflation.R")
  run_script("02_application/inflation/MZdiagnostic.R")

  financial_prepared_inputs <- c(
    "02_application/finrisk/data/emini_fcasts_vola.rds",
    "02_application/finrisk/data/emini_fcasts_VaR_1.rds",
    "02_application/finrisk/data/emini_fcasts_VaR_5.rds"
  )
  missing_financial_prepared_inputs <-
    financial_prepared_inputs[!file.exists(financial_prepared_inputs)]

  if (!length(missing_financial_prepared_inputs)) {
    run_script("02_application/finrisk/evaluation_vola.R")
    run_script("02_application/finrisk/evaluation_var.R")
    run_script("02_application/finrisk/timeseriesplot.R")
    run_script("02_application/finrisk/MZdiagnostic.R")
  } else {
    message(
      "\nPrepared financial data not found; financial outputs were skipped: ",
      paste(missing_financial_prepared_inputs, collapse = ", ")
    )
  }
}


message("\nSelected replication steps completed successfully.")
