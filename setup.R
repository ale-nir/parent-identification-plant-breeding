# =============================================================================
# setup.R
#
# Install (if necessary) and load all required packages.
# =============================================================================

required_packages <- c(
  "here",
  "dplyr",
  "stringr",
  "ggplot2",
  "ggpubr",
  "rstatix",
  "dunn.test",
  "Rcpp",
  "readr",
  "ggsignif"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {

  message("Installing missing packages...")
  install.packages(missing_packages, dependencies = TRUE)

}

invisible(
  lapply(required_packages, library, character.only = TRUE)
)


Rcpp::sourceCpp(here("src","jaccard_cpp.cpp"))
Rcpp::sourceCpp(here("src","lod_cpp.cpp"))

source(here("R","data_io.R"))
source(here("R","pedigree_processing.R"))

datasets <- load_project_data()

data <- datasets$data
bin_matrix <- datasets$bin_matrix

source(here("R","parentage_utils.R"))
source(here("R","jaccard_analysis.R"))
source(here("R","lod_analysis.R"))
