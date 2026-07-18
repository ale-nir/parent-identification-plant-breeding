# =============================================================================
# lod_analysis.R
#
# LOD score-based parentage analysis
#
# This file contains the functions required to identify putative parent–offspring
# relationships using pairwise LOD scores. The implemented workflow identifies
# the highest-ranking candidate parents for each evaluated variety and compares
# the retrieved candidates with documented pedigree information.
#
# The file provides functions to:
#   - identify the most likely related varieties based on LOD scores,
#   - evaluate the retrieval of documented parents,
#   - summarize the performance of the LOD-based parentage analysis.
#
# Common helper functions shared with other parentage analyses are implemented
# in parentage_utils.R.
# =============================================================================
# Identify the highest-ranking related varieties based on LOD scores.
#
# Retrieves the documented children and putative siblings of a selected variety
# using pedigree information. Siblings are defined as varieties sharing at
# least one documented parent, excluding the selected variety itself and its
# direct parents.
#
# Parameters
# ----------
# selected_variety : character
#     Name of the variety for which related genotypes should be identified.
#
# data : data.frame
#     Metadata table containing variety names and documented pedigree
#     information.
#
# Returns
# -------
# A list containing:
#
#   related
#       Union of children and siblings.
#
#   children
#       Varieties having the selected variety recorded as a parent.
#
#   siblings
#       Varieties sharing at least one documented parent with the selected
#       variety.
# =============================================================================

identify_related_varieties <- function(selected_variety,
                                       data){

  # Standardize names and pedigree information
  data <-
    data %>%
    dplyr::mutate(

      Name = tolower(Name),

      Available.Pedigree =
        .clean_pedigree(Available.Pedigree)

    )

  selected_variety <- tolower(selected_variety)

  # Retrieve documented children
  children <-

    data %>%
    dplyr::filter(

      stringr::str_detect(
        Available.Pedigree,
        stringr::fixed(selected_variety)
      )

    ) %>%
    dplyr::pull(Name)

  # Retrieve documented parents
  selected_variety_data <-

    data %>%
    dplyr::filter(Name == selected_variety)

  documented_parents <-

    if(nrow(selected_variety_data) > 0){

      as.character(selected_variety_data$Available.Pedigree)

    } else{

      ""

    }

  siblings <- character(0)

  if(length(documented_parents) > 0 &&
     all(!is.na(documented_parents)) &&
     all(nchar(documented_parents) > 0) &&
     all(documented_parents != "unknown")){

    parent_names <- .extract_parents(documented_parents)
	
	parent_names <-

	parent_names[
		!(tolower(parent_names) %in%
		c("", "na", "n/a", "unknown", "?", "null"))
	]

    for(parent in parent_names){

      siblings <-

        union(

          siblings,

          data$Name[
            stringr::str_detect(
              data$Available.Pedigree,
              stringr::fixed(parent)
            )
          ]

        )

    }

    siblings <- setdiff(siblings, selected_variety)
    siblings <- setdiff(siblings, parent_names)

  }

  related <- union(children, siblings)

  list(

    related  = toupper(related),
    children = toupper(children),
    siblings = toupper(siblings)

  )

}

# =============================================================================
# Compute summary statistics for the LOD-based parentage analysis.
#
# Computes the overall parent retrieval success, retrieval rates for each
# documented parent and descriptive statistics of the corresponding LOD
# scores from the detailed parentage analysis results.
#
# Parameters
# ----------
# parentage_results       : output list generated for each evaluated variety
# parent_scores_filtered  : dataframe containing successfully retrieved
#                           documented parent LOD scores
#
# Returns
# -------
# A list containing:
#   - success_rate_global
#   - success_rates_by_parent
#   - score_stats
# =============================================================================

.compute_summary_statistics_lod <- function(parentage_results,
                                            parent_scores_filtered){

  parent_stats <- data.frame(

    variety = names(parentage_results),

    parent1_in_top = sapply(
      parentage_results,
      function(x) x$parent1$in_top
    ),

    parent2_in_top = sapply(
      parentage_results,
      function(x) x$parent2$in_top
    ),

    stringsAsFactors = FALSE

  )

  parent_stats$any_parent_in_top <-
    parent_stats$parent1_in_top |
    parent_stats$parent2_in_top

  score_stats <- list(

    parent1 = list(

      n =
        sum(
          parent_scores_filtered$parent_type == "parent1"
        ),

      min =
        min(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent1"
          ],
          na.rm = TRUE
        ),

      max =
        max(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent1"
          ],
          na.rm = TRUE
        ),

      mean =
        mean(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent1"
          ],
          na.rm = TRUE
        ),

      median =
        median(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent1"
          ],
          na.rm = TRUE
        ),

      q1 =
        quantile(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent1"
          ],
          .25,
          na.rm = TRUE
        ),

      q3 =
        quantile(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent1"
          ],
          .75,
          na.rm = TRUE
        ),

      sd =
        sd(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent1"
          ],
          na.rm = TRUE
        )

    ),

    parent2 = list(

      n =
        sum(
          parent_scores_filtered$parent_type == "parent2"
        ),

      min =
        min(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent2"
          ],
          na.rm = TRUE
        ),

      max =
        max(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent2"
          ],
          na.rm = TRUE
        ),

      mean =
        mean(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent2"
          ],
          na.rm = TRUE
        ),

      median =
        median(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent2"
          ],
          na.rm = TRUE
        ),

      q1 =
        quantile(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent2"
          ],
          .25,
          na.rm = TRUE
        ),

      q3 =
        quantile(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent2"
          ],
          .75,
          na.rm = TRUE
        ),

      sd =
        sd(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent2"
          ],
          na.rm = TRUE
        )

    )

  )

  list(

    success_rate_global =
      mean(
        parent_stats$any_parent_in_top,
        na.rm = TRUE
      ),

    success_rates_by_parent = c(

      parent1 =
        mean(
          parent_stats$parent1_in_top,
          na.rm = TRUE
        ),

      parent2 =
        mean(
          parent_stats$parent2_in_top,
          na.rm = TRUE
        )

    ),

    score_stats = score_stats

  )

}

# =============================================================================
# Evaluate parent retrieval using LOD scores.
#
# Compares the highest-ranking LOD candidates with documented pedigree
# information to determine whether the recorded parents are successfully
# retrieved. Candidate parents are filtered according to the predefined
# exclusion criteria before ranking. Detailed retrieval results, parent
# score information and summary statistics describing the performance of
# the LOD-based analysis are returned.
#
# Parameters
# ----------
# data                     : dataframe containing variety metadata
# lod_matrix               : pairwise LOD score matrix
# varieties_with_parents   : dataframe of varieties with documented pedigree
# top_n                    : number of highest-ranking candidates to evaluate
#
# Returns
# -------
# A list containing:
#   - detailed_results        : retrieval results for each evaluated variety
#   - all_parent_scores       : LOD scores for all documented parents
#   - filtered_parent_scores  : documented parents successfully retrieved
#   - success_rate_global     : proportion of varieties with at least one
#                               documented parent retrieved
#   - success_rates_by_parent : retrieval rate for parent1 and parent2
#   - score_stats             : summary statistics of retrieved parent LOD
#                               scores
# =============================================================================

analyze_parentage_lod <- function(data,
                                  lod_matrix,
                                  varieties_with_parents,
                                  top_n){

  ## -------------------------------------------------------------------------
  ## Initialisation
  ## -------------------------------------------------------------------------

  parentage_results <- vector(
    "list",
    length(varieties_with_parents$Name)
  )

  names(parentage_results) <- varieties_with_parents$Name

  data_names_lower <- tolower(data$Name)

  name_to_index <- setNames(
    seq_len(nrow(data)),
    data_names_lower
  )

  parent_data <- .prepare_parent_data(
    data,
    varieties_with_parents
  )

  parent_score_records <- vector(
    "list",
    nrow(parent_data) * 2
  )

  record_index <- 1

  ## -------------------------------------------------------------------------
  ## Parentage analysis
  ## -------------------------------------------------------------------------

  for(i in seq_len(nrow(parent_data))){

    current_parent_info <- parent_data[i, ]

    variety_name <- current_parent_info$Name
    variety_index <- name_to_index[[tolower(variety_name)]]

    variety_year <- current_parent_info$Year

    lod_scores <- lod_matrix[variety_index, ]

    lod_df <- data.frame(
      Genotype = data$Name,
      Score_LOD = lod_scores,
      Year = data$Year,
      stringsAsFactors = FALSE
    )

    lod_df <- lod_df[
      lod_df$Genotype != variety_name,
    ]

    lod_df <- lod_df[
      lod_df$Year <= variety_year,
    ]

    related <- identify_related_varieties(
      variety_name,
      data
    )$related

    lod_df <- lod_df[
      !lod_df$Genotype %in% related,
    ]

    lod_df <- lod_df[
      order(-lod_df$Score_LOD),
    ]

    top_scores <- head(
      lod_df,
      top_n
    )

    top_candidates <- unname(
      name_to_index[
        tolower(top_scores$Genotype)
      ]
    )

    for(parent_type in c("parent1", "parent2")){

      parent_name <- current_parent_info[[parent_type]]

      if(is.na(parent_name))
        next

      parent_idx <- unname(
        name_to_index[tolower(parent_name)]
      )

      if(is.na(parent_idx))
        next

      parent_rank <- match(
        parent_idx,
        top_candidates
      )

      parent_score_records[[record_index]] <- data.frame(
        variety = variety_name,
        parent_type = parent_type,
        parent_name = parent_name,
        in_top = !is.na(parent_rank),
        rank = parent_rank,
        score = lod_scores[parent_idx],
        stringsAsFactors = FALSE
      )

      record_index <- record_index + 1

    }

    parentage_results[[variety_name]] <- list(

      variety = variety_name,

      all_scores = lod_df,

      top_scores = top_scores,

      excluded = c(
        related,
        data$Name[data$Year > variety_year]
      ),

      parent1 = .get_parent_info(
        current_parent_info$parent1,
        top_candidates,
        lod_scores,
        name_to_index
      ),

      parent2 = .get_parent_info(
        current_parent_info$parent2,
        top_candidates,
        lod_scores,
        name_to_index
      )

    )

  }

  ## -------------------------------------------------------------------------
  ## Parent score table
  ## -------------------------------------------------------------------------

  parent_scores <- do.call(
    rbind,
    parent_score_records[seq_len(record_index - 1)]
  )

  parent_scores_filtered <-
    parent_scores[
      !is.na(parent_scores$rank),
    ]

  ## -------------------------------------------------------------------------
  ## Summary statistics
  ## -------------------------------------------------------------------------

  summary_statistics <- .compute_summary_statistics_lod(
    parentage_results,
    parent_scores_filtered
  )

  ## -------------------------------------------------------------------------
  ## Return results
  ## -------------------------------------------------------------------------

  list(

    detailed_results = parentage_results,

    all_parent_scores = parent_scores,

    filtered_parent_scores = parent_scores_filtered,

    success_rate_global =
      summary_statistics$success_rate_global,

    success_rates_by_parent =
      summary_statistics$success_rates_by_parent,

    score_stats =
      summary_statistics$score_stats

  )

}