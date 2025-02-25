#!/usr/bin/env Rscript

options(crayon.enabled = TRUE)

# Set environment variables
Sys.setenv(LOGNAME = Sys.info()[["user"]])
Sys.setenv("_R_CHECK_FORCE_SUGGESTS_" = "false")
Sys.setenv("_R_CHECK_CRAN_INCOMING_" = "false")

# Define the check directory
check_dir_path <- file.path(getwd(), "check")

# Print equivalent output
cat("LOGNAME=", Sys.getenv("LOGNAME"), "\n")
cat("check-dir-path=", check_dir_path, "\n")

# Ensure required packages are installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("rcmdcheck", quietly = TRUE)) {
  install.packages("rcmdcheck")
}

# Delete NAMESPACE and man/ to force Roxygen to regenerate them
namespace_file <- "NAMESPACE"
man_dir <- "man"

if (file.exists(namespace_file)) {
  file.remove(namespace_file)
  cat("Deleted existing NAMESPACE file.\n")
}

if (dir.exists(man_dir)) {
  unlink(man_dir, recursive = TRUE)
  cat("Deleted existing man/ directory.\n")
}

# Rebuild documentation
cat("Running devtools::document()...\n")
devtools::document()

# Run package check
cat("Running rcmdcheck::rcmdcheck()...\n")
rcmdcheck::rcmdcheck(
  args = c("--no-manual", "--as-cran"),
  build_args = c("--no-manual", "--compact-vignettes=gs+qpdf"),
  error_on = "warning",
  check_dir = check_dir_path
)

BiocCheck::BiocCheck('new-package'=TRUE) 
