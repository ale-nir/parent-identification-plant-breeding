# ==============================================================================
# combined_analysis.R
#
# Combined analysis of Jaccard similarity and LOD parentage inference
#
# This file provides functions to compare the results obtained from Jaccard
# similarity and LOD score analyses. It identifies concordant and discordant
# parent recoveries for documented parent–offspring relationships and summarizes
# the agreement between both metrics. The functions support either a Top-k
# selection strategy or a Jaccard score interval, allowing a direct comparison
# of parent retrieval performance under different selection criteria.
# ==============================================================================
# Compute summary statistics for the combined Jaccard–LOD analysis
#
# Calculates the number and percentage of documented parent recoveries obtained
# by both methods, by Jaccard only, by LOD only, or by neither method. The
# function also reports overall concordant and discordant classifications based
# on the combined analysis summary table.
#
# Args:
#   summary_table: Data frame returned by analyze_parentage_combined()
#                  containing the recovery status for each documented parent.
#
# Returns:
#   A data frame summarizing absolute counts and percentages for each recovery
#   category.
# ==============================================================================

.compute_combined_statistics <- function(summary_table){

  n <- nrow(summary_table)

  both <- sum(summary_table$Recovered_By == "Both")
  jaccard_only <- sum(summary_table$Recovered_By == "Jaccard")
  lod_only <- sum(summary_table$Recovered_By == "LOD")
  neither <- sum(summary_table$Recovered_By == "None")

  concordant <- both + neither
  discordant <- jaccard_only + lod_only

  stats <- data.frame(
    Category = c(
      "Both",
      "Jaccard only",
      "LOD only",
      "Neither",
      "Concordant",
      "Discordant"
    ),
    Count = c(
      both,
      jaccard_only,
      lod_only,
      neither,
      concordant,
      discordant
    ),
    Percentage = round(
      100 * c(
        both,
        jaccard_only,
        lod_only,
        neither,
        concordant,
        discordant
      ) / n,
      2
    ),
    stringsAsFactors = FALSE
  )

  stats
}

# ------------------------------------------------------------------------------
# Compare parent retrieval between Jaccard similarity and LOD scores
#
# Performs a pairwise comparison of documented parent recoveries obtained from
# Jaccard similarity and LOD analyses. For each known parent–offspring
# relationship, the function records whether the parent is recovered by each
# method according to either a Top-k ranking or a Jaccard score interval. It
# produces a detailed comparison table together with summary statistics
# describing the agreement between both approaches.
#
# Args:
#   jaccard_results: Result object returned by analyze_parentage_jaccard().
#   lod_results: Result object returned by analyze_parentage_lod().
#   jaccard_selection: Selection strategy for Jaccard results. Either "top"
#                      (Top-k ranking) or "range" (score interval).
#   top_k: Maximum accepted rank when using the Top-k strategy.
#   jaccard_min: Lower bound of the accepted Jaccard score interval when
#                jaccard_selection = "range".
#   jaccard_max: Upper bound of the accepted Jaccard score interval when
#                jaccard_selection = "range".
#
# Returns:
#   A list containing:
#     - summary_table: Detailed comparison for each documented parent,
#       including Jaccard and LOD ranks, scores, recovery status, and
#       concordance.
#     - global_stats: Summary of absolute counts and percentages for each
#       recovery category (Both, Jaccard only, LOD only, Neither, Concordant,
#       Discordant).
# ------------------------------------------------------------------------------

analyze_parentage_combined <- function(
    jaccard_results,
    lod_results,
    jaccard_selection = c("top", "range"),
    top_k = 10,
    jaccard_min = NULL,
    jaccard_max = NULL
){

  jaccard_selection <- match.arg(jaccard_selection)

  if(jaccard_selection == "range"){

    if(is.null(jaccard_min) || is.null(jaccard_max))
      stop(
        "jaccard_min and jaccard_max must be provided ",
        "when jaccard_selection = 'range'."
      )

  }

  ## ----------------------------------------------------------
  ## Retrieve detailed results
  ## ----------------------------------------------------------

  jaccard_details <- jaccard_results$detailed_results
  lod_details <- lod_results$detailed_results

  ## ----------------------------------------------------------
  ## Keep only varieties available in both analyses
  ## ----------------------------------------------------------

  common_varieties <-
    intersect(
      names(jaccard_details),
      names(lod_details)
    )

  out <- vector("list", length(common_varieties) * 2)

  k <- 1

  ## ----------------------------------------------------------
  ## Loop over varieties
  ## ----------------------------------------------------------

  for(var_name in common_varieties){

    jac <- jaccard_details[[var_name]]
    lod <- lod_details[[var_name]]

    for(parent in c("parent1", "parent2")){

      jac_parent <- jac[[parent]]
      lod_parent <- lod[[parent]]

      if(is.null(jac_parent) ||
         is.null(lod_parent) ||
         is.na(jac_parent$name))
        next

      ## ----------------------------------------------------------
      ## Jaccard selection
      ## ----------------------------------------------------------

      if(jaccard_selection == "top"){

        selected_j <-
          !is.na(jac_parent$rank) &&
          jac_parent$rank <= top_k

      }else{

        selected_j <-
          !is.na(jac_parent$score) &&
          jac_parent$score >= jaccard_min &&
          jac_parent$score <= jaccard_max

      }

      ## ----------------------------------------------------------
      ## LOD selection
      ## ----------------------------------------------------------

      selected_lod <-
        !is.na(lod_parent$rank) &&
        lod_parent$rank <= top_k

      ## ----------------------------------------------------------
      ## Store results
      ## ----------------------------------------------------------

      out[[k]] <-
        data.frame(
          Variety = var_name,
          Parent = jac_parent$name,
          Rank_J = jac_parent$rank,
          Score_J = jac_parent$score,
          Rank_LOD = lod_parent$rank,
          Score_LOD = lod_parent$score,
          Selected_J = selected_j,
          Selected_LOD = selected_lod,
          stringsAsFactors = FALSE
        )

      k <- k + 1

    }

  }

  summary_table <-
    do.call(
      rbind,
      out[seq_len(k - 1)]
    )

  ## ----------------------------------------------------------
  ## Recovery categories
  ## ----------------------------------------------------------

  summary_table$Recovered_By <-
    ifelse(
      summary_table$Selected_J & summary_table$Selected_LOD,
      "Both",
      ifelse(
        summary_table$Selected_J,
        "Jaccard",
        ifelse(
          summary_table$Selected_LOD,
          "LOD",
          "None"
        )
      )
    )

  summary_table$Discordant <-
    xor(
      summary_table$Selected_J,
      summary_table$Selected_LOD
    )
	
  global_stats <- .compute_combined_statistics(summary_table)

  ## ----------------------------------------------------------
  ## Summary statistics
  ## ----------------------------------------------------------
	return(
  list(
    summary_table = summary_table,
	global_stats = global_stats
   
  ))

}