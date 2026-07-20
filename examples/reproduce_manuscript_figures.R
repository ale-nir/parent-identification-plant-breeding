# ==============================================================================
# Reproduce the figures presented in the manuscript
#
# This script illustrates the workflow used to generate the datasets required
# for the publication figures.
#
# Several intermediate objects produced during this workflow (e.g. parentage
# analysis results and relationship score datasets) are also the source of the
# summary tables reported in the manuscript.
#
# The individual functions used here are documented separately in the
# documentation files located in the docs/ directory.
# ==============================================================================

source("setup.R")

# ==============================================================================
# Figure 1
# Jaccard similarity distribution for documented parent-offspring pairs.
#
# The parentage analysis results generated in this section are also used to
# compute the parent retrieval statistics reported in the manuscript.
# ==============================================================================

ped <- create_pedigree_subsets(data)

varieties_both_parents <- ped$biparental_subset
varieties_any_parent   <- ped$uniparental_subset

J <- compute_jaccard_matrix_cpp(bin_matrix)

# parentage_both and parentage_any contain:
# - ranked candidate parents
# - parent retrieval statistics
# - documented parent scores (used for Figure 1)

parentage_both <- analyze_parentage_jaccard(
    data,
    J,
    varieties_both_parents,
    top_k = 10
)

# Example of parental global retrieval statistics
# parentage_both$success_rate_global

# Example of detailed statistics
# parentage_both$score_stats

parentage_any <- analyze_parentage_jaccard(
    data,
    J,
    varieties_any_parent,
    top_k = 10
)

parent_pairs_both <- parentage_both$all_parent_scores
parent_pairs_any  <- parentage_any$all_parent_scores

bi_pairs <- data.frame(
    jaccard.score = parent_pairs_both$score,
    type = "Bi-parental"
)

uni_pairs <- data.frame(
    jaccard.score = parent_pairs_any$score,
    type = "Single-parent"
)

plot_data <- dplyr::bind_rows(
    bi_pairs,
    uni_pairs
)

p <- plot_parental_zone_histogram(plot_data)

save_publication_figure(
  plot = p,
  filename = "jaccard_distribution_plot"
)

# ==============================================================================
# Figure 2
# Structural comparison of relationship groups.
#
# The relationship score dataset generated here is also used for the
# descriptive statistics and non-parametric tests reported in the manuscript.
# ==============================================================================

# relationship_both and relationship_any contain:
# - parent, sibling and random relationship scores
# - descriptive statistics
# - the data used for the structural comparison (Figure 2)

relationship_both  <- generate_jaccard_relationship_scores(
    data,
    J,
    varieties_both_parents
)

relationship_any <- generate_jaccard_relationship_scores(
    data,
    J,
    varieties_any_parent
)

relationship_both$Group <- factor(
    relationship_both$Group,
    levels = c("Parent", "Sibling", "Random")
)

relationship_any$Group <- factor(
    relationship_any$Group,
    levels = c("Parent", "Sibling", "Random")
)

relationship_both$Subset <- "Biparental"
relationship_any$Subset <- "Single-parent"

# Retrieve the information on any relationship pair 
# head(relationship_both)

plot_data <- dplyr::bind_rows(
    relationship_both,
    relationship_any
)

p <- plot_relationship_boxplot(plot_data)

save_publication_figure(
  plot = p,
  filename = "jaccard_structural_analysis"
)