# ==============================================================================
# plotting_utils.R
#
# Utility functions for generating publication-quality figures from the outputs
# of the parentage analyses.
#
# This file contains reusable ggplot2 functions for visualizing Jaccard and LOD
# score distributions, parent retrieval performance and relationship-based
# comparisons. The functions are designed to reproduce the figures presented in
# the accompanying manuscript while remaining reusable for other datasets.
#
# The functions only generate graphical objects and do not perform the analyses
# themselves. Statistical analyses should first be carried out using the
# functions provided in jaccard_analysis.R, lod_analysis.R and
# combined_analysis.R.
# ==============================================================================
# Plot Jaccard similarity distributions highlighting the parental zone
#
# Generate a histogram and kernel density estimate of Jaccard similarity scores
# for one or more pedigree subsets. A user-defined parental similarity zone is
# highlighted to facilitate visual assessment of parent-offspring similarity
# distributions.
#
# Parameters
# ----------
# data : data.frame
#     Data frame containing at least the columns:
#         jaccard.score : Jaccard similarity scores.
#         type          : Pedigree subset labels.
#
# parental_zone : numeric
#     Numeric vector of length two defining the lower and upper bounds of the
#     highlighted parental similarity zone.
#
# binwidth : numeric
#     Histogram bin width.
#
# Returns
# -------
# A ggplot object.
# ==============================================================================
plot_parental_zone_histogram <- function(data,
                                         parental_zone = c(0.61, 0.69),
                                         binwidth = 0.02) {

  ggplot(data, aes(x = jaccard.score)) +

    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = binwidth,
      fill = "grey85",
      color = "black"
    ) +

    geom_density(
      color = "black",
      linewidth = 0.3,
      adjust = 1.5
    ) +

    annotate(
      "rect",
      xmin = parental_zone[1],
      xmax = parental_zone[2],
      ymin = 0,
      ymax = Inf,
      alpha = 0.15,
      fill = "steelblue"
    ) +

    geom_vline(
      xintercept = parental_zone,
      linetype = "dashed",
      linewidth = 0.3
    ) +

    scale_x_continuous(
      breaks = seq(0.35, 0.80, by = 0.05)
    ) +

    facet_wrap(
      ~type,
      scales = "fixed",
      ncol = 1
    ) +

    labs(
      x = "Jaccard similarity",
      y = "Density"
    ) +

    theme_minimal(base_size = 13)

}

# ==============================================================================
# Plot Jaccard similarity distributions by relationship group
#
# Generate boxplots comparing Jaccard similarity scores among documented
# parents, siblings and unrelated varieties for one or more pedigree subsets.
#
# Parameters
# ----------
# data : data.frame
#     Data frame containing at least the following columns:
#         Group   : relationship category.
#         Jaccard : Jaccard similarity score.
#         Subset  : pedigree subset label.
#
# Returns
# -------
# A ggplot object.
# ==============================================================================
plot_relationship_boxplot <- function(data) {

  ggplot(
    data,
    aes(
      x = Group,
      y = Jaccard,
      fill = Group
    )
  ) +

    geom_boxplot(
      width = 0.6,
      outlier.shape = NA,
      alpha = 0.85
    ) +

    geom_jitter(
      width = 0.15,
      alpha = 0.15,
      size = 0.6
    ) +

    scale_fill_manual(
      values = c(
        "#5DA5A4",
        "#C7A900",
        "#E67E5F"
      )
    ) +

    facet_wrap(
      ~Subset,
      scales = "fixed"
    ) +

    labs(
      x = "Relationship category",
      y = "Jaccard similarity"
    ) +

    coord_cartesian(
      ylim = c(0.30, 0.88)
    ) +

    theme_minimal(base_size = 14) +

    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.spacing = unit(1.5, "lines"),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.7
      )
    )
}

# ==============================================================================
# Save a ggplot figure in publication quality
#
# Save a ggplot object in JPEG or TIFF format using publication-quality
# settings.
#
# Parameters
# ----------
# plot : ggplot
#     ggplot object to save.
#
# filename : character
#     Output filename without extension.
#
# format : character
#     Output image format. Either "jpeg" or "tiff".
#
# width : numeric
#     Figure width.
#
# height : numeric
#     Figure height.
#
# units : character
#     Size units passed to ggsave().
#
# dpi : numeric
#     Image resolution.
#
# Returns
# -------
# The saved file path (invisibly).
# ==============================================================================
save_publication_figure <- function(plot,
                                    filename,
                                    format = c("jpeg", "tiff"),
                                    width = 29,
                                    height = 23,
                                    units = "cm",
                                    dpi = 600) {

  format <- match.arg(format)

  if (format == "jpeg") {

    ggsave(
      filename = paste0(filename, ".jpeg"),
      plot = plot,
      width = width,
      height = height,
      units = units,
      dpi = dpi,
      device = "jpeg",
      bg = "white"
    )

  } else {

    ggsave(
      filename = paste0(filename, ".tiff"),
      plot = plot,
      width = width,
      height = height,
      units = units,
      dpi = dpi,
      device = "tiff",
      compression = "lzw",
      bg = "white"
    )

  }

  invisible(
    paste0(filename, ".", ifelse(format == "jpeg", "jpeg", "tiff"))
  )
}