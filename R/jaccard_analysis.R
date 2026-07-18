# =============================================================================
# jaccard_analysis.R
#
# Functions for ranking candidate parents using Jaccard similarity scores.
#
# The Jaccard similarity matrix is assumed to be precomputed (Rcpp).
# This file performs:
#
#   - parent retrieval
#   - ranking analysis
#   - Top-k evaluation
#   - summary statistics
#
# =============================================================================
# Compute summary statistics for the parent retrieval analysis.
#
# Parameters
# ----------
# parentage_results       : output list generated for each evaluated variety
# parent_scores_filtered  : dataframe containing retrieved parent scores
#
# Returns
# -------
# A list containing:
#   - success_rate_global
#   - success_rates_by_parent
#   - score_stats
# =============================================================================

.compute_summary_statistics_jaccard <- function(parentage_results,
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

  list(

    success_rate_global =
      mean(parent_stats$any_parent_in_top,
           na.rm = TRUE),

    success_rates_by_parent = c(

      parent1 =
        mean(parent_stats$parent1_in_top,
             na.rm = TRUE),

      parent2 =
        mean(parent_stats$parent2_in_top,
             na.rm = TRUE)

    ),

    score_stats = list(

      parent1 =
        summary(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent1"
          ]
        ),

      parent2 =
        summary(
          parent_scores_filtered$score[
            parent_scores_filtered$parent_type == "parent2"
          ]
        )

    )

  )
}

# =============================================================================
# analyze_parentage_jaccard
#
# Evaluates the ability of Jaccard similarity scores to retrieve documented
# parents among the Top-k ranked candidate genotypes.
#
# Parameters
# ----------
# data : metadata dataframe.
#
# jaccard_matrix : symmetric Jaccard similarity matrix.
#
# varieties_with_parents : dataframe containing the varieties included in the
#                          parentage analysis.
#
# top_n : number of highest-ranked candidates retained.
#
# Returns
# -------
# A list containing:
#
#   detailed_results
#       Detailed ranking for each evaluated variety.
#
#   all_parent_scores
#       Jaccard scores of all documented parents.
#
#   filtered_parent_scores
#       Parent scores restricted to retrieved parents.
#
#   success_rate_global
#       Fraction of varieties for which at least one documented parent was
#       retrieved among the Top-k candidates.
#
#   success_rates_by_parent
#       Retrieval rate for parent 1 and parent 2 separately.
#
#   score_stats
#       Summary statistics of retrieved parent scores.
# =============================================================================

analyze_parentage_jaccard <- function(data,
                                      jaccard_matrix,
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

    variety_name <- parent_data$Name[i]
    variety_index <- name_to_index[[tolower(variety_name)]]

    jaccard_scores <- jaccard_matrix[variety_index, ]

    ranking <- order(
      jaccard_scores,
      decreasing = TRUE
    )

    ranking <- ranking[ranking != variety_index]

    top_candidates <- ranking[seq_len(top_n)]

    top_scores <- data.frame(
      Genotype = data$Name[top_candidates],
      Score_Jaccard = jaccard_scores[top_candidates],
      stringsAsFactors = FALSE
    )

    all_scores <- data.frame(
      Genotype = data$Name[ranking],
      Score_Jaccard = jaccard_scores[ranking],
      stringsAsFactors = FALSE
    )

    current_parent_info <- parent_data[i, ]

    for(parent_type in c("parent1", "parent2")){

      parent_name <- current_parent_info[[parent_type]]

      if(is.na(parent_name))
        next

      parent_idx <- unname(name_to_index[parent_name])

      if(is.na(parent_idx))
        next

      parent_rank <- match(parent_idx, top_candidates)

      parent_score_records[[record_index]] <- data.frame(
        variety = variety_name,
        parent_type = parent_type,
        parent_name = parent_name,
        in_top = !is.na(parent_rank),
        rank = parent_rank,
        score = jaccard_scores[parent_idx],
        stringsAsFactors = FALSE
      )

      record_index <- record_index + 1

    }

    parentage_results[[variety_name]] <- list(

      variety = variety_name,

      all_scores = all_scores,

      top_scores = top_scores,

      parent1 = .get_parent_info(
        current_parent_info$parent1,
        top_candidates,
        jaccard_scores,
        name_to_index
      ),

      parent2 = .get_parent_info(
        current_parent_info$parent2,
        top_candidates,
        jaccard_scores,
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
    parent_scores[!is.na(parent_scores$rank), ]

  ## -------------------------------------------------------------------------
  ## Summary statistics
  ## -------------------------------------------------------------------------

  summary_statistics <- .compute_summary_statistics_jaccard(
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