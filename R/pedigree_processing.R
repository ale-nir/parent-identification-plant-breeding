# =============================================================================
# pedigree_processing.R
#
# Functions for processing pedigree information and creating the pedigree
# subsets used throughout the manuscript.
#
# The pedigree is expected to be stored in the column
# "Available.Pedigree" using the format:
#
# Parent1 x Parent2
#
# =============================================================================

.clean_pedigree <- function(x){

  x %>%
    tolower() %>%
    gsub("\\([^)]*\\)", "", .) %>%
    gsub("\\s+", " ", .) %>%
    stringr::str_trim()

}

.extract_parents <- function(pedigree){

  .clean_pedigree(pedigree) %>%
    strsplit(" x ") %>%
    .[[1]] %>%
    stringr::str_trim()

}

.classify_pedigree <- function(parents, all_names){

  p1_present <- length(parents) >= 1 &&
                !is.na(parents[1]) &&
                parents[1] %in% all_names

  p2_present <- length(parents) >= 2 &&
                !is.na(parents[2]) &&
                parents[2] %in% all_names

  if(p1_present && p2_present)
    return("biparental")

  if(p1_present || p2_present)
    return("uniparental")

  return("none")

}

create_pedigree_subsets <- function(data){

  all_names <- tolower(data$Name)

  pedigree_class <- character(nrow(data))

  for(i in seq_len(nrow(data))){

    parents <- .extract_parents(data$Available.Pedigree[i])

    pedigree_class[i] <- .classify_pedigree(
      parents,
      all_names
    )

  }

  list(

    biparental_subset =
      data[pedigree_class == "biparental", ],

    uniparental_subset =
      data[pedigree_class == "uniparental", ]

  )

}

