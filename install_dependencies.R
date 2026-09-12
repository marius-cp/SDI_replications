script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  setwd(dirname(normalizePath(sub("^--file=", "", script_arg))))
}

reference_r_version <- "4.2.1"
if (as.character(getRversion()) != reference_r_version) {
  stop(
    "Use R ", reference_r_version,
    " for the reference environment (current version: ", getRversion(), ").",
    call. = FALSE
  )
}

cran_repository <- "https://cloud.r-project.org"

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = cran_repository)
}

cran_packages <- c(
  car = "3.1.1",
  doParallel = "1.0.17",
  dplyr = "1.0.10",
  fGarch = "4033.92",
  foreach = "1.5.2",
  geomtextpath = "0.1.1",
  ggh4x = "0.2.8",
  ggnewscale = "0.4.10",
  ggplot2 = "3.5.2",
  ggtext = "0.1.2",
  gridExtra = "2.3",
  gtable = "0.3.6",
  kableExtra = "1.3.4",
  knitr = "1.40",
  lubridate = "1.8.0",
  midasr = "0.8",
  murphydiagram = "0.12.2",
  patchwork = "1.3.1",
  purrr = "0.3.5",
  quantreg = "5.94",
  readr = "2.1.3",
  readxl = "1.4.1",
  rugarch = "1.4.8",
  sandwich = "3.1.1",
  tibble = "3.1.8",
  tidyr = "1.2.1",
  tidyverse = "1.3.2",
  zoo = "1.8.12"
)

for (package_name in names(cran_packages)) {
  expected_version <- unname(cran_packages[[package_name]])
  installed_version <- if (requireNamespace(package_name, quietly = TRUE)) {
    as.character(utils::packageVersion(package_name))
  } else {
    NA_character_
  }

  if (is.na(installed_version) || installed_version != expected_version) {
    remotes::install_version(
      package = package_name,
      version = expected_version,
      repos = cran_repository,
      upgrade = "never"
    )
  }
}

expected_sdi_sha <- "0edf793ea183061d3d6eb53cbb60faaffb5b23e9"
installed_sdi_sha <- if (requireNamespace("SDI", quietly = TRUE)) {
  utils::packageDescription("SDI")$RemoteSha
} else {
  NULL
}

if (is.null(installed_sdi_sha) || installed_sdi_sha != expected_sdi_sha) {
  remotes::install_github(
    paste0("marius-cp/SDI@", expected_sdi_sha),
    upgrade = "never"
  )
}

source("00_functions/check_dependencies.R")
