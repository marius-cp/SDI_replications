required_r_version <- "4.2.1"

required_packages <- c(
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
  SDI = "0.0.1",
  tibble = "3.1.8",
  tidyr = "1.2.1",
  tidyverse = "1.3.2",
  zoo = "1.8.12"
)

notes <- character()
missing_packages <- character()

if (as.character(getRversion()) != required_r_version) {
  notes <- c(
    notes,
    sprintf("R %s is installed; the reference run used R %s.", getRversion(), required_r_version)
  )
}

for (package_name in names(required_packages)) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    missing_packages <- c(missing_packages, package_name)
    next
  }

  installed_version <- as.character(utils::packageVersion(package_name))
  expected_version <- unname(required_packages[[package_name]])
  if (installed_version != expected_version) {
    notes <- c(
      notes,
      sprintf(
        "Package '%s' is version %s; the reference run used %s.",
        package_name,
        installed_version,
        expected_version
      )
    )
  }
}

if (requireNamespace("SDI", quietly = TRUE)) {
  expected_sdi_sha <- "0edf793ea183061d3d6eb53cbb60faaffb5b23e9"
  installed_sdi_sha <- utils::packageDescription("SDI")$RemoteSha
  if (is.null(installed_sdi_sha) || installed_sdi_sha != expected_sdi_sha) {
    notes <- c(
      notes,
      sprintf(
        "SDI is not recorded at the reference commit %s (installed record: %s).",
        expected_sdi_sha,
        if (is.null(installed_sdi_sha)) "missing" else installed_sdi_sha
      )
    )
  }
}

if (!isTRUE(capabilities("cairo"))) {
  notes <- c(notes, "R was built without Cairo PDF support; PDF output may fail or differ.")
}

if (length(missing_packages)) {
  stop(
    "Required package(s) are not installed: ",
    paste(missing_packages, collapse = ", "),
    ". See REPRODUCIBILITY.md.",
    call. = FALSE
  )
}

if (length(notes)) {
  message(
    paste(
      c(
        "The software environment differs from the reference environment; the replication will continue:",
        paste0("- ", notes),
        "Reference versions are recorded in SESSION_INFO.txt."
      ),
      collapse = "\n"
    )
  )
} else {
  message("All installed versions match the reference environment (R ", required_r_version, ").")
}
