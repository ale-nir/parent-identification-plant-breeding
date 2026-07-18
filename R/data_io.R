# =============================================================================
# data_io.R
#
# Load project datasets
#
# Reads the metadata and binary genotype matrices from CSV files.
# The first column of the binary file is assumed to contain sample identifiers
# and is converted to row names.
#
# Parameters
# ----------
# data_file : path to metadata CSV file
# bin_file  : path to binary matrix CSV file
#
# Returns
# -------
# A list containing:
#   - data       : metadata dataframe
#   - bin_matrix : binary genotype matrix
# =============================================================================

load_project_data <- function(data_file = here("data", "data.csv"),
                              bin_file  = here("data", "bin.csv")) {

  data <- read.csv(
    data_file,
    sep = ";",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  bin <- read.csv(
    bin_file,
    sep = ";",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  data <- data[,-1]

  rownames(bin) <- bin[[1]]
  bin <- bin[, -1, drop = FALSE]

  bin_matrix <- as.matrix(bin)
  mode(bin_matrix) <- "numeric"

  return(list(
    data = data,
    bin_matrix = bin_matrix
  ))
}