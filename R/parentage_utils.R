# =============================================================================
# Prepare pedigree information for the parentage analysis.
#
# Extracts the documented parents for each evaluated variety and stores them
# in two dedicated columns (parent1 and parent2). Parent names are converted
# to lower case to facilitate subsequent matching with the genotype database.
#
# Parameters
# ----------
# data : metadata dataframe containing pedigree information.
#
# varieties_with_parents : dataframe containing the varieties included in the
#                          parentage analysis.
#
# Returns
# -------
# A dataframe containing one row per evaluated variety with the additional
# columns:
#
#   parent1
#       First documented parent.
#
#   parent2
#       Second documented parent (if available).
# =============================================================================

.prepare_parent_data <- function(data,
                                 varieties_with_parents){

  parent_data <-
    data %>%
    dplyr::filter(Name %in% varieties_with_parents$Name)
    

  parents <- lapply(parent_data$Available.Pedigree, .extract_parents)

  parent_data$parent1 <-
    sapply(
      parents,
      function(x)
        if(length(x) >= 1) x[1] else NA_character_
    )

  parent_data$parent2 <-
    sapply(
      parents,
      function(x)
        if(length(x) >= 2) x[2] else NA_character_
    )

  parent_data

}

# =============================================================================
# Retrieve information for a documented parent.
#
# Parameters
# ----------
# parent_name     : name of the documented parent
# top_candidates  : indices of the Top-k ranked candidates
# score_vector    : similarity scores for the query variety
# name_to_index   : named vector mapping genotype names to row indices
#
# Returns
# -------
# A list containing:
#   - name
#   - in_top
#   - rank
#   - score
# =============================================================================

.get_parent_info <- function(parent_name,
                             top_candidates,
                             score_vector,
                             name_to_index){

  if(is.na(parent_name)){
    return(list(
      name = NA_character_,
      in_top = FALSE,
      rank = NA_integer_,
      score = NA_real_
    ))
  }

  parent_idx <- unname(name_to_index[tolower(parent_name)])

  if(is.na(parent_idx)){
    return(list(
      name = parent_name,
      in_top = FALSE,
      rank = NA_integer_,
      score = NA_real_
    ))
  }

  parent_rank <- match(parent_idx, top_candidates)

  list(
    name = parent_name,
    in_top = !is.na(parent_rank),
    rank = parent_rank,
    score = score_vector[parent_idx]
  )

}
